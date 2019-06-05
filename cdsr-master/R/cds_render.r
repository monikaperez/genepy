#' Render a RMarkdown file and upload it to https://dev.cds.team/cds_projects
#'
#' @importFrom rmarkdown render
#' @importFrom RCurl postForm fileUpload
#' @importFrom jsonlite fromJSON
#' @param path_to_rmd The path to the RMarkdown file you would want to render and upload to cds_projects
#' @param folder Folder to store your file into. Can't be an already Git taken folder name
#' @param path_to_html The path to the HTML File you already generated on your owm. We will skip the rendering of the RMD, so please make sure both match
#' @param ... All parameteres you would want to pass to the rmarkdown render function
#' @return The full path of the HTML generated
#' @examples
#' cds_render("local_rmd.Rmd")
#' @export cds_render
cds_render <- function(path_to_rmd, folder=NULL, path_to_html=NULL, ...) {
    if(!is.null(path_to_html)) {
      print("### Using the html provided ###")
      if(!file.exists(path_to_html)){
        stop(paste("File", path_to_html, "does not exist. Are you sure about the path provided?"))
      }
    } else {
      print("### Rendering RMarkdown ###")
      path_to_html <- rmarkdown::render(path_to_rmd, ...)
    }


    # Upload to CDS_Projects
    print("### Uploading to CDS_Projects ###")
    #RCurl::postForm("https://dev.cds.team/cds_projects/api/upload", html=RCurl::fileUpload(path_to_html), rmd=RCurl::fileUpload(path_to_rmd))
    tryCatch(postReturn <- RCurl::postForm("https://dev.cds.team/cds_projects/api/upload",
                                  html=RCurl::fileUpload(path_to_html),
                                  rmd=RCurl::fileUpload(path_to_rmd),
                                  folder_name=folder),
             error = function(x){
                 message <- x$message
                 if(message == 'FORBIDDEN\r\n'){
                    stop(paste0(folder, ' is a git folder. Please use another name or add your html using git'))
                 }else{
                    stop(x)
                 }
             })

    r_postReturn <- jsonlite::fromJSON(postReturn)
    print(paste0("### Link to the project is ", r_postReturn$projectUrl, " ###"))
    print(paste0("### Link to your file is ", r_postReturn$fileUrl, " ###"))
    return(path_to_html)
}
