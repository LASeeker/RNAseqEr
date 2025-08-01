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