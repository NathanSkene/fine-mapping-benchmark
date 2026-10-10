# Release integrity and custody are separate from scientific blindness.

.challenge_sha256 <- function(path) {
  exe <- Sys.which("openssl")
  if (!nzchar(exe)) stop("openssl is required for SHA-256 release integrity")
  out <- system2(exe, c("dgst", "-sha256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  hit <- regmatches(out, regexpr("[a-fA-F0-9]{64}$", out))
  if (length(hit) != 1L || nchar(hit) != 64L) stop("Could not hash artifact")
  tolower(hit)
}

#' Write separately stored public inputs and private assessor truth
#' @param panel Output of generate_challenge().
#' @param public_dir New directory containing inputs and public metadata only.
#' @param private_dir New directory for truth and generator parameters. It must
#'   be outside public_dir; a live blind round needs independent access control.
#' @param release_id Dataset release identifier distinct from specification/run.
#' @return Public and private manifest paths with SHA-256 file commitments.
#' @export
write_challenge_release <- function(panel, public_dir, private_dir,
                                     release_id = "development-v1") {
  for (s in list(public_dir, private_dir, release_id))
    if (!is.character(s) || length(s) != 1L || is.na(s) || !nzchar(s)) stop("Invalid release paths/ID")
  # Expand relative paths before checking custody separation.
  absolute <- function(path) {
    path <- path.expand(path)
    if (!grepl("^(/|[A-Za-z]:)", path)) path <- file.path(getwd(), path)
    if (grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", path)) stop("Paths must not contain '..'")
    sub("/+$", "", normalizePath(path, mustWork = FALSE))
  }
  pub <- absolute(public_dir); priv <- absolute(private_dir)
  if (pub == priv || startsWith(priv, paste0(pub, "/")) || startsWith(pub, paste0(priv, "/")))
    stop("Public and private release directories must be separate")
  if (file.exists(pub) || file.exists(priv)) stop("Refusing to overwrite a release")
  .challenge_validate_panel(panel)
  for (x in panel$inputs) validate_challenge_input(x)
  dir.create(pub, recursive = TRUE); dir.create(priv, recursive = TRUE, mode = "0700")
  # Verify post-creation resolved paths too (including any symlink parents).
  pub <- normalizePath(pub); priv <- normalizePath(priv)
  if (startsWith(priv, paste0(pub, "/")) || startsWith(pub, paste0(priv, "/")) || pub == priv)
    stop("Resolved release paths overlap")
  saveRDS(panel$inputs, file.path(pub, "inputs.rds"), version = 2)
  metadata <- list(specification = panel$config$specification,
                   release_id = release_id, n_tasks = length(panel$inputs),
                   split = panel$config$split,
                   contract = "public-input-v1", truth_available = FALSE)
  saveRDS(metadata, file.path(pub, "metadata.rds"), version = 2)
  saveRDS(panel$truth, file.path(priv, "truth.rds"), version = 2)
  saveRDS(panel$config, file.path(priv, "generator.rds"), version = 2)
  manifest <- function(directory, files) {
    d <- data.frame(file = files, sha256 = vapply(file.path(directory, files), .challenge_sha256, ""))
    utils::write.csv(d, file.path(directory, "manifest.csv"), row.names = FALSE)
    file.path(directory, "manifest.csv")
  }
  list(public_manifest = manifest(pub, c("inputs.rds", "metadata.rds")),
       private_manifest = manifest(priv, c("truth.rds", "generator.rds")))
}

#' Verify a frozen release and load predictor inputs
#' @param public_dir Directory created by write_challenge_release().
#' @return List with public inputs and metadata, after integrity/schema checks.
#' @export
read_challenge_release <- function(public_dir) {
  m <- utils::read.csv(file.path(public_dir, "manifest.csv"), stringsAsFactors = FALSE)
  if (!identical(names(m), c("file", "sha256")) || anyDuplicated(m$file) ||
      !setequal(m$file, c("inputs.rds", "metadata.rds"))) stop("Invalid manifest")
  found <- vapply(file.path(public_dir, m$file), .challenge_sha256, "")
  if (!identical(unname(found), m$sha256)) stop("Artifact hash mismatch")
  inputs <- readRDS(file.path(public_dir, "inputs.rds"))
  metadata <- readRDS(file.path(public_dir, "metadata.rds"))
  for (x in inputs) validate_challenge_input(x)
  if (metadata$n_tasks != length(inputs) || is.null(names(inputs)) || anyDuplicated(names(inputs)) ||
      !identical(names(inputs), unname(vapply(inputs, `[[`, "", "task_id")))) stop("Invalid release identities")
  list(inputs = inputs, metadata = metadata)
}

#' Adapt an existing region wrapper to the public challenge contract
#' @param wrapper Function taking region_geno, region_pheno and method arguments.
#' @param method_args Prespecified wrapper arguments; never truth or tuning labels.
#' @return Predictor function. It uses summary statistics only and explicitly
#'   rejects multi-cohort or variable-n inputs needing specialised adapters.
#' @export
challenge_region_adapter <- function(wrapper, method_args = list()) {
  if (!is.function(wrapper) || !is.list(method_args)) stop("Invalid wrapper arguments")
  function(input) {
    validate_challenge_input(input)
    if (!is.null(input$cohorts) || !is.null(input$n_per_variant))
      stop("This region adapter supports a single constant-n cohort")
    p <- length(input$z)
    geno <- list(n = input$n, p = p, LD = input$LD, LD_ref = input$LD,
                 maf = rep(.25, p), variant_ids = input$variant_id)
    pheno <- list(n = input$n, z = input$z, beta_hat = input$beta_hat,
                  se = input$se, LD = input$LD, annotations_matrix = input$annotations)
    result <- do.call(wrapper, c(list(region_geno = geno, region_pheno = pheno), method_args))
    if (!is.null(result$error)) stop(result$error)
    native <- result$credible_sets %||% list()
    sets <- lapply(native, function(s) {
      if (!is.numeric(s) || any(!is.finite(s)) || any(s != floor(s)) || any(s < 1 | s > p))
        stop("Wrapper returned invalid credible-set indices")
      input$variant_id[s]
    })
    list(task_id = input$task_id, variant_id = input$variant_id,
         pip = result$pip, credible_sets = sets, status = "ok",
         implementation = result$method, oracle = FALSE,
         converged = result$additional$converged %||% NA)
  }
}
