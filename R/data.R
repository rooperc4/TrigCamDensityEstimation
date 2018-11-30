#' bocaccio

#' This is a data set containing data collected at the Footprint Bank (Channel Islands, California) 
#' in October 2017 using stationary Triggered Cameras. The data are zero-filled counts of bocaccio
#' rockfish by frame (captured at 30 second intevals). Images were analyzed using SEBASTES software
#' (Williams et al. 2016). These collections were made during the Untrawlable Habitat Strategic 
#' Initiative. The methods are documented in Williams et al. 2018.

#' @format A data frame with 770130 rows and 11 variables:
#' \describe{
#'   \item{deployment_ID}{unique identifier for the deployment}
#'   \item{frame_number}{individual frame number}
#'   \item{range_bin}{bin in 0.5 m increments from 1 m from the camera to 6 m from camera}
#'   \item{unit}{TrigCam unit number}
#'   \item{site}{deployment site}
#'   \item{frame_time}{time image was taken (UTC)}
#'   \item{primary_habitat_type}{substrate class}
#'   \item{secondary_habitat_type}{substrate class}
#'   \item{volume}{volume of the range_bin estimated from bottom detection}
#'   \item{class}{fish species (NA's for where there was zero density)}
#'   \item{density}{density in range bin and frame}
#'   ...
#' }
#' @source {UHSI project}
"bocaccio"




