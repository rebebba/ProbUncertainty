#' Star Trek The Original Series - Survival of Starfleet personnel
#'
#' A dataset containing the survival status (dead or alive) for each of the 430 crew 
#' members on-board the Enterprise at the end of season 3 and their Starfleet uniform 
#' color. The crew included 43 officers (gold uniforms), 151 scientists (blue uniforms),
#' and 236 members of the operations division (red uniforms).
#'
#' @format A data frame with 430 rows and 2 variables. Each row is a unique crew member.
#' \describe{
#'   \item{uniform}{the color of the Starfleet uniform (gold, blue, or red) worn by crew member i}
#'   \item{status}{the survival status at the end of season 3 (dead or alive) of crew member i}
#' }
#' @source <https://www.ex-astris-scientia.org/database/redshirt_deaths.htm>
#' @source <https://memory-gamma.fandom.com/wiki/USS_Enterprise_(NCC-1701)>
#' 
#' @usage data("startrek")
"startrek"