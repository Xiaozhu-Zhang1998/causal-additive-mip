# Shared PNS, ordering, pruning, bootstrap and RESIT wrappers.
## CAM ========

# Step 1: CAM preliminary neighborhood selection ====

cam_pns <- function(
    X,
    at_least_selected = 0.1,
    max_neighbors = 10L,
    output = FALSE
) {
  
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (n < 2L || p < 2L) {
    stop("X must contain at least two rows and two columns.")
  }
  
  if (anyNA(X) || any(!is.finite(X))) {
    stop("CAM-PNS requires finite data without missing values.")
  }
  
  pars <- list(
    atLeastThatMuchSelected = at_least_selected,
    atMostThatManyNeighbors = as.integer(max_neighbors)
  )
  
  sel_gam_boost <- getFromNamespace("selGamBoost", "CAM")
  
  ## Column j = candidate parents of node j
  sel_mat <- sapply(seq_len(p), function(j) {
    as.logical(
      sel_gam_boost(
        X = X,
        pars = pars,
        output = output,
        k = j
      )
    )
  })
  
  diag(sel_mat) <- FALSE
  
  rownames(sel_mat) <- colnames(X)
  colnames(sel_mat) <- colnames(X)
  
  ## Directed edge list
  directed_edges <- which(sel_mat, arr.ind = TRUE)
  
  edge_table <- data.frame(
    parent = colnames(X)[directed_edges[, 1]],
    child = colnames(X)[directed_edges[, 2]],
    parent_index = directed_edges[, 1],
    child_index = directed_edges[, 2],
    row.names = NULL
  )
  
  ## Undirected OR skeleton
  skeleton_or <- sel_mat | t(sel_mat)
  diag(skeleton_or) <- FALSE
  
  list(
    candidate_parents = sel_mat,
    candidate_skeleton = skeleton_or,
    edge_table = edge_table,
    parameters = pars
  )
}


fixed_superstructure_selection <- function(
    X,
    pars,
    output = FALSE,
    k
) {
  selMat <- as.matrix(pars$selMat)
  p <- ncol(X)
  
  if (!identical(dim(selMat), c(p, p))) {
    stop("selMat must be a p x p matrix.")
  }
  
  selected <- as.logical(selMat[, k])
  selected[k] <- FALSE
  
  if (output) {
    cat(
      "Candidate parents for node", k, ":",
      which(selected), "\n"
    )
  }
  
  selected
}

# Step 2: Fitting CAM ====

cam_estimate_order <- function(
    X,
    selMat,
    score_name = "SEMGAM",
    num_basis_fcts = 10L,
    max_num_parents = NULL,
    num_cores = 1L,
    output = FALSE
) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  p <- ncol(X)
  n <- nrow(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  
  selMat <- as.matrix(selMat)
  
  if (!identical(dim(selMat), c(p, p))) {
    stop("selMat must have dimension p x p.")
  }
  
  selMat <- selMat != 0
  diag(selMat) <- FALSE
  
  if (is.null(max_num_parents)) {
    max_num_parents <- min(
      p - 1L,
      floor(n / 20)
    )
  }
  
  fit <- CAM::CAM(
    X = X,
    scoreName = score_name,
    parsScore = list(
      numBasisFcts = as.integer(num_basis_fcts)
    ),
    numCores = as.integer(num_cores),
    maxNumParents = as.integer(max_num_parents),
    output = output,
    
    variableSel = TRUE,
    variableSelMethod = fixed_superstructure_selection,
    variableSelMethodPars = list(
      selMat = selMat
    ),
    
    # Turn off Step 3
    pruning = FALSE,
    
    intervData = FALSE,
    intervMat = NULL
  )
  
  unpruned_adj <- as.matrix(fit$Adj)
  
  dimnames(unpruned_adj) <- list(
    parent = colnames(X),
    child = colnames(X)
  )
  
  order_index <- topological_order_from_adj(unpruned_adj)
  
  list(
    order = order_index,
    order_names = colnames(X)[order_index],
    unpruned_adjacency = unpruned_adj,
    score = fit$Score,
    raw_fit = fit
  )
}

topological_order_from_adj <- function(adj) {
  adj <- as.matrix(adj)
  
  if (nrow(adj) != ncol(adj)) {
    stop("adj must be square.")
  }
  
  p <- nrow(adj)
  A <- adj != 0
  diag(A) <- FALSE
  
  indegree <- colSums(A)
  available <- which(indegree == 0L)
  
  order <- integer(0)
  
  while (length(available) > 0L) {
    # Deterministic tie-breaking
    node <- min(available)
    available <- setdiff(available, node)
    
    order <- c(order, node)
    
    children <- which(A[node, ])
    
    if (length(children) > 0L) {
      A[node, children] <- FALSE
      indegree[children] <- indegree[children] - 1L
      
      newly_available <- children[indegree[children] == 0L]
      available <- sort(unique(c(available, newly_available)))
    }
  }
  
  if (length(order) != p) {
    stop("The supplied adjacency matrix contains a directed cycle.")
  }
  
  order
}



# Step 3: Variable selection from ordering ====

gam_prune_from_order <- function(
    X,
    order,
    selMat = NULL,
    cutoff = 0.001,
    num_basis_fcts = 10L,
    output = FALSE
) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  p <- ncol(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  
  order <- as.integer(order)
  
  if (
    length(order) != p ||
    !setequal(order, seq_len(p))
  ) {
    stop("order must be a permutation of 1, ..., p.")
  }
  
  if (is.null(selMat)) {
    selMat <- matrix(TRUE, p, p)
    diag(selMat) <- FALSE
  } else {
    selMat <- as.matrix(selMat)
    
    if (!identical(dim(selMat), c(p, p))) {
      stop("selMat must have dimension p x p.")
    }
    
    selMat <- selMat != 0
    diag(selMat) <- FALSE
  }
  
  sel_gam <- getFromNamespace("selGam", "CAM")
  
  final_adj <- matrix(
    0L,
    nrow = p,
    ncol = p,
    dimnames = list(
      parent = colnames(X),
      child = colnames(X)
    )
  )
  
  for (position in seq_len(p)) {
    child <- order[position]
    
    if (position == 1L) {
      next
    }
    
    predecessors <- order[seq_len(position - 1L)]
    
    # Enforce the common directed superstructure
    candidate_parents <- predecessors[
      selMat[predecessors, child]
    ]
    
    if (length(candidate_parents) == 0L) {
      next
    }
    
    if (output) {
      cat(
        "Child:", colnames(X)[child],
        "\nCandidate parents:",
        colnames(X)[candidate_parents],
        "\n"
      )
    }
    
    # selGam expects the response to be the kth/final column
    local_data <- cbind(
      X[, candidate_parents, drop = FALSE],
      X[, child]
    )
    
    selected_local <- sel_gam(
      X = local_data,
      pars = list(
        cutOffPVal = cutoff,
        numBasisFcts = as.integer(num_basis_fcts)
      ),
      output = output,
      k = ncol(local_data)
    )
    
    selected_parent_positions <- which(
      selected_local[seq_along(candidate_parents)]
    )
    
    if (length(selected_parent_positions) > 0L) {
      selected_parents <-
        candidate_parents[selected_parent_positions]
      
      final_adj[selected_parents, child] <- 1L
    }
  }
  
  final_adj
}


## NPVAR ========
fit_cam_once = function(X) {
  cam_order_result <- suppressWarnings(
    cam_estimate_order(
      X = X,
      selMat = selMat,
      output = FALSE
    )
  )
  cam_order <- cam_order_result$order
  
  cam_adj <- gam_prune_from_order(
    X = X,
    order = cam_order,
    selMat = selMat,
    cutoff = 0.001,
    num_basis_fcts = 10,
    output = FALSE
  )
  return(cam_adj)
}


# fit npvar ----
fit_npvar_once = function(X) {
  npvar_order <- NPVAR(
    X,
    run.np = FALSE,
    run.mgcv = TRUE,
    layer.select = FALSE
  )
  
  npvar_adj <- gam_prune_from_order(
    X = X,
    order = npvar_order,
    selMat = selMat,
    cutoff = 0.001,
    num_basis_fcts = 10,
    output = FALSE
  )
  return(npvar_adj)
}


# bootstrap ----
bootstrap_cam_npvar <- function(
    X,
    fit_cam_once,
    fit_npvar_once,
    B = 100,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  cam_counts <- matrix(0, p, p)
  npvar_counts <- matrix(0, p, p)
  
  for (b in seq_len(B)) {
    index <- sample.int(n, size = n, replace = TRUE)
    X_boot <- X[index, , drop = FALSE]
    
    cam_adj <- fit_cam_once(X_boot)
    cam_counts <- cam_counts + (cam_adj != 0)
    
    npvar_adj <- fit_npvar_once(X_boot)
    npvar_counts <- npvar_counts + (npvar_adj != 0)
    
    cat("finished", b, '\n')
  }
  
  cam_prob_adj <- cam_counts / B
  diag(cam_prob_adj) <- 0
  
  npvar_prob_adj <- npvar_counts / B
  diag(npvar_prob_adj) <- 0
  
  return(list(
    cam_prob_adj=cam_prob_adj,
    npvar_prob_adj=npvar_prob_adj
  ))
}


## EqVAR ========
resit_estimate_order <- function(
    M,
    alpha = 0.05,
    model = train_gam,
    parsModel = list(numBasisFcts = 10),
    indtest = indtestHsic,
    parsIndtest = list(method = "ExactFastTrace"),
    output = FALSE
) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  
  if (anyNA(M) || any(!is.finite(M))) {
    stop("M must contain only finite, nonmissing values.")
  }
  
  n <- nrow(M)
  p <- ncol(M)
  
  if (n < 2L || p < 2L) {
    stop("M must contain at least two observations and two variables.")
  }
  
  if (is.null(colnames(M))) {
    colnames(M) <- paste0("X", seq_len(p))
  }
  
  remaining <- seq_len(p)
  sink_order <- integer(p)
  
  # Store diagnostics from each elimination step
  pvalue_history <- vector("list", p - 1L)
  
  for (step in seq_len(p - 1L)) {
    candidate_pvalues <- numeric(length(remaining))
    
    for (k in seq_along(remaining)) {
      candidate_sink <- remaining[k]
      regressors <- remaining[-k]
      
      if (output) {
        cat(
          "Step", step,
          ": fitting", colnames(M)[candidate_sink],
          "using",
          paste(colnames(M)[regressors], collapse = ", "),
          "\n"
        )
      }
      
      test_result <- fit_and_test_independence(
        x = M[, candidate_sink],
        y = M[, regressors, drop = FALSE],
        z = M[, regressors, drop = FALSE],
        alpha = alpha,
        model = model,
        parsModel = parsModel,
        indtest = indtest,
        parsIndtest = parsIndtest
      )
      
      candidate_pvalues[k] <- test_result$p.value
      
      if (output) {
        cat(
          "  HSIC p-value:",
          format(test_result$p.value, digits = 5),
          "\n"
        )
      }
    }
    
    # Original ICML code uses:


    #
    # This is equivalent to selecting the largest p-value.
    best_position <- which.max(candidate_pvalues)
    selected_sink <- remaining[best_position]
    
    sink_order[step] <- selected_sink
    
    pvalue_history[[step]] <- data.frame(
      step = step,
      candidate_index = remaining,
      candidate = colnames(M)[remaining],
      p_value = candidate_pvalues,
      selected = remaining == selected_sink,
      row.names = NULL
    )
    
    if (output) {
      cat(
        "Selected sink:",
        colnames(M)[selected_sink],
        "with p-value",
        format(candidate_pvalues[best_position], digits = 5),
        "\n\n"
      )
    }
    
    remaining <- remaining[-best_position]
  }
  
  # The last remaining node is the source-most node.
  sink_order[p] <- remaining[1L]
  
  # ICML stores sink -> source, so reverse it.
  topological_order <- rev(sink_order)
  
  list(
    order = topological_order,
    order_names = colnames(M)[topological_order],
    sink_order = sink_order,
    sink_order_names = colnames(M)[sink_order],
    pvalue_history = pvalue_history,
    alpha = alpha,
    model = deparse(substitute(model)),
    parsModel = parsModel,
    indtest = deparse(substitute(indtest)),
    parsIndtest = parsIndtest
  )
}
