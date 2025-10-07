#' @title Variable Selection (Forward or Backward) for models of \code{BSW()}
#'
#' @description
#' Performs forward or backward variable selection based on Wald test p-values for models estimated using \code{bsw()}.
#' In each step, a new model is fitted using \code{bsw()}, and variables are added or removed based on the significance level defined by \code{alpha}.
#'
#' @usage variable_selection_bsw(model, selection = c("backward", "forward"), alpha = 0.157,
#'                               print_models = FALSE, maxit = NULL, conswitch = NULL)
#'
#' @param model A model object from \code{bsw()} with full data and formula.
#' @param selection Character string, either \code{"backward"} or \code{"forward"}. Determines the direction of model selection. If not specified, backward elimination is performed by default.
#' @param alpha P-value threshold for variable inclusion (forward) or exclusion (backward). Defaults to 0.157, as recommended by Heinze, G., Wallisch, C., & Dunkler, D. (2018).
#' @param print_models Logical; whether to print each model during selection. Defaults to FALSE.
#' @param maxit Maximum number of iterations in the bsw() algorithm. If NULL, defaults to 200L or value from original model call.
#' @param conswitch Specifies how the constraint matrix is constructed:
#' \describe{
#'   \item{1 (default)}{Generates all possible combinations of minimum and maximum values for the predictors (excluding the intercept), resulting in \eqn{2^{m-1}} constraints.
#'   This formulation constrains model predictions within the observed data range, making it suitable for both risk factor identification and prediction (prognosis).}
#'   \item{0}{Uses the raw design matrix \code{x} as the constraint matrix, resulting in \eqn{n} constraints.
#'   This is primarily suitable for identifying risk factors, but not for prediction tasks, as predictions are not bounded to realistic ranges.}
#' }
#' @references Heinze, G., Wallisch, C., & Dunkler, D. (2018). Variable selection – A review and recommendations for the practicing statistician. Biometrical Journal, 60(3), 431–449.
#' @importFrom stats model.frame model.matrix as.formula terms
#'
#' @return An object of class \code{"bsw_selection"}, which is a list containing:
#' \describe{
#'   \item{final_model}{An object of class \code{bsw} representing the final model selected through the variable selection process.}
#'   \item{model_list}{A list of intermediate \code{bsw} model objects fitted during each step of the selection.}
#'   \item{skipped_models}{A named list of models that failed to converge and were skipped during the selection. Each entry includes the attempted formula.}
#'   \item{final_formula}{The final model formula used in the last step.}
#'   \item{EPV}{Estimated events-per-variable (EPV) of the final model, used as a diagnostic for model stability.}
#'   \item{warnings}{Optional warning messages about convergence issues or model stability (e.g., low EPV or skipped variables).}
#' }
#'
#' @author Julius Johannes Weise, Thomas Wolf, Stefan Wagenpfeil
#' @examples
#' set.seed(123)
#' x1 <- rnorm(500, 50, 10)
#' x2 <- rnorm(500, 30, 5)
#' x3 <- rnorm(500, 40, 8)
#' x4 <- rnorm(500, 60, 12)
#' logit <- (-4 + x1 * 0.04 + x3 * 0.04)
#' p <- 1 / (1 + exp(-logit))
#' y <- rbinom(500, 1, p)
#' df <- data.frame(y, x1, x2, x3, x4)
#' fit <- bsw(formula = y ~ x1 + x2 + x3 + x4, data = df)
#' result <- variable_selection_bsw(fit, selection = "forward", alpha = 0.1)
#' print(result)
#'
#' @export

variable_selection_bsw <- function(model, selection = c("backward", "forward"),
                                   alpha = 0.157, print_models = FALSE,
                                   maxit = NULL, conswitch = NULL) {
  # Input validation
  selection <- match.arg(selection)
  

  if (!inherits(model, "bsw")) {
    stop("\"model\" must be an object of class \"bsw\" as returned by bsw()")
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("\"alpha\" must be a single numeric value between 0 and 1 (exclusive)")
  }

  if (!is.logical(print_models) || length(print_models) != 1L) {
    stop("\"print_models\" must be TRUE or FALSE")
  }

  if (!is.null(maxit) && (!is.numeric(maxit) || length(maxit) != 1L || maxit <= 0)) {
    stop("\"maxit\" must be a single positive number or NULL")
  }

  if (!is.null(conswitch) && !conswitch %in% c(0, 1)) {
    stop("\"conswitch\" must be 0 or 1")
  }

  data <- model@data
  full_formula <- model@formula

  terms_obj <- terms(full_formula, data = data)
  response <- as.character(attr(terms_obj, "variables")[[2]])
  all_predictors <- attr(terms_obj, "term.labels")

  call_list <- as.list(model@call)

  # Warnungsliste vorbereiten
  warn_text <- character()

  if (is.null(maxit)) {
    maxit <- if (!is.null(call_list$maxit)) eval(call_list$maxit) else 200L
  }
  
  if (is.null(conswitch)) {
    conswitch <- if (!is.null(call_list$conswitch)) eval(call_list$conswitch) else 1
  }

  model_list <- list()
  model_step <- 1

  skipped_models <- list()

  #Backward Selection
  if (selection == "backward") {

    predictors <- all_predictors

    while (length(predictors) > 0) {
      current_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
      current_model <- bsw(current_formula, data, maxit = maxit, conswitch = conswitch)

      model_list[[length(model_list) + 1]] <- current_model

      if (!current_model@converged) {
        warning("Model did not converge. Stopping selection.")
        break
      }

      if (print_models) {
        cat("\n--- Model step:", model_step, "---\n")
        try(summary(current_model), silent = TRUE)
      }

      model_step <- model_step + 1

      coefs <- coef(current_model)
      se <- sqrt(diag(solve(hess(coefs, current_model@y, current_model@x))))
      z <- coefs / se
      p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      p <- p[-1]

      worst_var <- names(p)[which.max(p)]

      if (p[worst_var] > alpha) {
        if (length(predictors) > 1) {
          predictors <- setdiff(predictors, worst_var)
        } else {
          warn_text <- c(warn_text, "No variable remained below the p-value threshold.")
          break
        }
      }
      else {
        break
      }
    }

    final_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  }

  # Forward Selection
  else if (selection == "forward") {

    selected_predictors <- c()
    remaining_predictors <- all_predictors

    repeat {
      p_values <- c()
      candidate_models <- list()

      for (var in remaining_predictors) {
        test_formula <- as.formula(paste(response, "~", paste(c(selected_predictors, var), collapse = " + ")))

        test_model <- tryCatch({
          bsw(test_formula, data, maxit = maxit, conswitch = conswitch)
        }, error = function(e) {
          return(NULL)
        })

        coefs <- coef(test_model)
        se <- tryCatch({
          sqrt(diag(solve(hess(coefs, test_model@y, test_model@x))))
        }, error = function(e) {
          return(NULL)
        })

        z <- coefs / se
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

        p_values[var] <- p[var]
        candidate_models[[var]] <- test_model

        if (is.null(test_model) || is.null(se) || !test_model@converged) {
          skipped_models[[paste0("Model step ", model_step, " variable ", var)]] <- list(
            formula = test_formula
          )
          next
        }
      }

      if (length(p_values) == 0) {
        if (length(remaining_predictors) == 0) {
          warn_text <- c(warn_text, "All variables have already been included.")
        } else {
          warn_text <- c(warn_text, "No further variables could be added due to failed convergence of all candidate models.")
        }
        break
      }


      best_var <- names(p_values)[which.min(p_values)]

      if (p_values[best_var] <= alpha) {
        selected_predictors <- c(selected_predictors, best_var)
        remaining_predictors <- setdiff(remaining_predictors, best_var)
        best_model <- candidate_models[[best_var]]
        model_list[[length(model_list) + 1]] <- best_model

        if (print_models) {
          cat("\n--- Model step:", model_step, "---\n")
          try(summary(best_model), silent = TRUE)
        }
        model_step <- model_step + 1
      } else {
        break
      }
    }

    # Warning if there are no predictors above the threshold
    if (length(selected_predictors) == 0) {
      warning("No variable met the significance criterion or all models failed to converge")
      return(invisible(NULL))
    }

    final_formula <- as.formula(paste(response, "~", paste(selected_predictors, collapse = " + ")))
  }

  final_model <- bsw(final_formula, data, maxit = maxit, conswitch = conswitch)


  # Stability Investigations
  y <- data[[response]]
  num_events <- sum(y == 1)
  num_predictors <- length(coef(final_model)) - 1
  epv <- num_events / max(1, num_predictors)


  # Warnung bei niedriger EPV
  if (epv < 25) {
    warn_text <- c(warn_text, "EPV < 25: Model may be unstable.")
  }

  # Warnung bei übersprungenen Modellen
  if (length(skipped_models) > 0) {
    skipped_names <- names(skipped_models)
    warn_text <- c(
      warn_text,
      paste0(length(skipped_models), " model(s) were skipped due to non-convergence."))
  }

  # Ergebnisobjekt vorbereiten
  out <- list(
    final_model = final_model,
    model_list = model_list,
    skipped_models = skipped_models,
    final_formula = final_formula,
    EPV = epv,
    warnings = if (length(warn_text) > 0) paste(warn_text, collapse = "\n") else NULL
  )

  class(out) <- "bsw_selection"
  invisible(out)

}
