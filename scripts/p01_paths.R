p01_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }
  getwd()
}

p01_project_root <- function() {
  env_root <- Sys.getenv("P01_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = FALSE))
  }
  normalizePath(file.path(p01_script_dir(), ".."), mustWork = FALSE)
}

# Mirrors p01_paths.source_data_root(): <ResearchD>/Data unless P01_SOURCE_DATA_ROOT is set.
p01_source_data_root <- function() {
  env_root <- Sys.getenv("P01_SOURCE_DATA_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = FALSE))
  }
  normalizePath(file.path(p01_project_root(), "..", "..", "Data"), mustWork = FALSE)
}
