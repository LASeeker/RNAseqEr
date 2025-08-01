#' Load environment variables from .env file using dotenv package
#' @description
#' Loads environment variables from a .env file using the dotenv package
#'
#' @param env_file Path to .env file (default: ".env" in current directory)
#'
#' @return Logical indicating if environment variables were loaded successfully
#' @export
#'
#' @examples
#' load_env_vars()
#' load_env_vars(".env")
#'
load_env_vars <- function(env_file = ".env") {
  # Check if dotenv package is available
  if(!requireNamespace("dotenv", quietly = TRUE)) {
    warning("dotenv package not installed. Installing it now...")
    install.packages("dotenv")
    if(!requireNamespace("dotenv", quietly = TRUE)) {
      stop("Failed to install dotenv package. Please install it manually: install.packages('dotenv')")
    }
  }
  
  # Load the dotenv package
  library(dotenv)
  
  if(file.exists(env_file)) {
    # Load environment variables using dotenv
    load_dot_env(file = env_file)
    
    message("✅ Environment variables loaded from ", env_file, " using dotenv package")
    return(TRUE)
  } else {
    warning("Environment file ", env_file, " not found")
    return(FALSE)
  }
}

#' Get API key from environment
#' @description
#' Retrieves API key from environment variables with fallback options
#'
#' @param key_name Name of the environment variable (e.g., "API_KEY_GPT4")
#' @param provider LLM provider for context-specific key names
#'
#' @return API key as character string, or NULL if not found
#' @export
#'
#' @examples
#' get_api_key("API_KEY_GPT4")
#' get_api_key(provider = "openai")
#' get_api_key(provider = "anthropic")
#'
get_api_key <- function(key_name = NULL, provider = NULL) {
  
  # If key_name is provided, use it directly
  if(!is.null(key_name)) {
    api_key <- Sys.getenv(key_name)
    if(api_key != "") {
      return(api_key)
    }
  }
  
  # If provider is specified, try provider-specific keys
  if(!is.null(provider)) {
    provider_keys <- list(
      "openai" = c("API_KEY_GPT4", "OPENAI_API_KEY", "GPT4_API_KEY"),
      "anthropic" = c("ANTHROPIC_API_KEY", "CLAUDE_API_KEY"),
      "local" = c("LOCAL_LLM_API_KEY", "LOCAL_API_KEY")
    )
    
    if(provider %in% names(provider_keys)) {
      for(key in provider_keys[[provider]]) {
        api_key <- Sys.getenv(key)
        if(api_key != "") {
          return(api_key)
        }
      }
    }
  }
  
  # Try common API key names
  common_keys <- c("API_KEY_GPT4", "OPENAI_API_KEY", "ANTHROPIC_API_KEY", 
                   "API_KEY", "LLM_API_KEY")
  
  for(key in common_keys) {
    api_key <- Sys.getenv(key)
    if(api_key != "") {
      return(api_key)
    }
  }
  
  return(NULL)
}

#' Check API key availability
#' @description
#' Checks if API keys are available for different LLM providers
#'
#' @param providers Vector of providers to check (default: all supported)
#'
#' @return List with availability status for each provider
#' @export
#'
#' @examples
#' check_api_keys()
#' check_api_keys(c("openai", "anthropic"))
#'
check_api_keys <- function(providers = c("openai", "anthropic", "local")) {
  results <- list()
  
  for(provider in providers) {
    api_key <- get_api_key(provider = provider)
    results[[provider]] <- !is.null(api_key)
  }
  
  return(results)
}

#' Setup environment for RNAseqEr
#' @description
#' Loads environment variables and checks API key availability
#'
#' @param env_file Path to .env file
#' @param auto_load Whether to automatically load environment variables
#'
#' @return List with setup status and available providers
#' @export
#'
#' @examples
#' setup_rnaseqer_env()
#' setup_rnaseqer_env(auto_load = FALSE)
#'
setup_rnaseqer_env <- function(env_file = ".env", auto_load = TRUE) {
  
  # Load environment variables if requested
  env_loaded <- FALSE
  if(auto_load) {
    env_loaded <- load_env_vars(env_file)
  }
  
  # Check API key availability
  api_status <- check_api_keys()
  
  # Provide helpful messages
  if(env_loaded) {
    message("✅ Environment variables loaded successfully using dotenv")
  } else {
    message("⚠️  Environment variables not loaded. Check if .env file exists.")
  }
  
  available_providers <- names(api_status)[unlist(api_status)]
  if(length(available_providers) > 0) {
    message("✅ Available LLM providers: ", paste(available_providers, collapse = ", "))
  } else {
    message("⚠️  No API keys found. LLM integration will not be available.")
  }
  
  return(list(
    env_loaded = env_loaded,
    api_status = api_status,
    available_providers = available_providers
  ))
}

#' Set API key manually for testing
#' @description
#' Manually sets the API key for testing the RNAseqEr workflow
#'
#' @param api_key Your OpenAI API key
#' @param provider LLM provider (default: "openai")
#'
#' @return Logical indicating if the key was set successfully
#' @export
#'
#' @examples
#' set_api_key_manual("your-api-key-here")
#'
set_api_key_manual <- function(api_key, provider = "openai") {
  if(provider == "openai") {
    Sys.setenv("API_KEY_GPT4" = api_key)
    Sys.setenv("OPENAI_API_KEY" = api_key)
  } else if(provider == "anthropic") {
    Sys.setenv("ANTHROPIC_API_KEY" = api_key)
    Sys.setenv("CLAUDE_API_KEY" = api_key)
  }
  
  # Verify it was set
  test_key <- get_api_key(provider = provider)
  if(!is.null(test_key)) {
    cat("✅ API key set successfully for", provider, "\n")
    cat("Key length:", nchar(test_key), "\n")
    return(TRUE)
  } else {
    cat("❌ Failed to set API key for", provider, "\n")
    return(FALSE)
  }
}

#' Test API key setup
#' @description
#' Tests if the API key is working properly
#'
#' @param api_key Your API key (optional, will use environment if not provided)
#' @param provider LLM provider (default: "openai")
#'
#' @return List with test results
#' @export
#'
#' @examples
#' test_api_key_setup()
#' test_api_key_setup("your-api-key-here")
#'
test_api_key_setup <- function(api_key = NULL, provider = "openai") {
  
  # Set API key if provided
  if(!is.null(api_key)) {
    set_api_key_manual(api_key, provider)
  }
  
  # Test environment setup
  env_status <- setup_rnaseqer_env(auto_load = FALSE)
  
  # Test API key retrieval
  retrieved_key <- get_api_key(provider = provider)
  
  # Test all providers
  all_providers <- check_api_keys()
  
  results <- list(
    env_loaded = env_status$env_loaded,
    api_key_found = !is.null(retrieved_key),
    api_key_length = if(!is.null(retrieved_key)) nchar(retrieved_key) else 0,
    available_providers = names(all_providers)[unlist(all_providers)],
    all_providers = all_providers
  )
  
  # Print results
  cat("=== API Key Test Results ===\n")
  cat("Environment loaded:", results$env_loaded, "\n")
  cat("API key found:", results$api_key_found, "\n")
  if(results$api_key_found) {
    cat("API key length:", results$api_key_length, "\n")
  }
  cat("Available providers:", paste(results$available_providers, collapse = ", "), "\n")
  
  return(results)
} 