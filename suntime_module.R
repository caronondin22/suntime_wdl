#suntime model wdl module script
#Author: Caroline Nondin (cnondin@fredhutch.org)
#date created: 6/4/26
#date last modified: 6/8/26

#Libraries
library(dplyr)
library(lubridate)
library(optparse)

option_list <- list(
  make_option("--tile_path",          type = "character", help = "Path to input tile .rds file"),
  make_option("--border_points_path", type = "character", help = "Path to border points .csv file"),
  make_option("--year",               type = "integer",   help = "Year for solar calculations (e.g. 2022)"),
  make_option("--matched_prefix",     type = "character",  default = "matched_", help = "Filename prefix for matched results (default: 'matched_')")
)

opt <- parse_args(OptionParser(option_list = option_list))

tile_path          <- opt$tile_path
border_points_path <- opt$border_points_path
year               <- opt$year
matched_prefix     <- opt$matched_prefix

# Derive tile name from input filename for flexible naming
tile_name <- tools::file_path_sans_ext(basename(tile_path))

#grabbing datasets we need
#border points
border_points <- read.csv(border_points_path)
required_cols <- c("timezone", "sunrise_avg_tz", "sunset_avg_tz")
missing_cols <- setdiff(required_cols, names(border_points))
if (length(missing_cols) > 0) {
  stop("Border points file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

#Upload tile
cat("\n Uploading tile:", tile_name, "\n")
tzpoints_time <- readRDS(tile_path)

## calculating sunrise and sunset time
cat("\n Calculating sunrise/sunset time \n")

## convert 43200 seconds to 12:00:00 o'clock just for displaying
## for calculation below, need to use the seconds version of time "time_clean")
tzpoints_time$time <- seconds_to_period(43200)
tzpoints_time$loopnum <- 365
tzpoints_time$date_start <- format(as.Date(paste0(year, "-01-01")), "%m/%d/%Y")
tzpoints_time$date_end <- format(as.Date(paste0(year, "-12-31")), "%m/%d/%Y")
tzpoints_time$year <- year
tzpoints_time$date_start_clean <- as.numeric(as.Date(paste0(year, "-01-01")) - as.Date("1960-01-01"))
tzpoints_time$date_end_clean <- as.numeric(as.Date(paste0(year, "-12-31")) - as.Date("1960-01-01"))
tzpoints_time$time_clean <- as.numeric(lubridate::hms(tzpoints_time$time))
tzpoints_time$sunrise_sum <- 0 #sunrise sum for the year
tzpoints_time$sunset_sum <- 0 #sunset sum for the year
tzpoints_time$valid_days <- 0 #counter for nb of valid days (has sunset/ sunrise values) per year 

#########################################################
## Step 2: calculating NOAA Solar Calculator variables ##
#########################################################

for (i in 1:unique(tzpoints_time$loopnum)) {

  tzpoints_time$pi <- pi

  current_date <- tzpoints_time$date_start_clean + (i - 1)

  ### julian day
  tzpoints_time$julianday <- round(current_date + 2415018.5 + tzpoints_time$time_clean - tzpoints_time$timezone/24,4)

  ### julian century
  tzpoints_time$juliancentury <- (tzpoints_time$julianday-2451545)/36525

  ### geom mean long sun (deg)
  tzpoints_time$geom_long <- (280.46646 + tzpoints_time$juliancentury*(36000.76983 + tzpoints_time$juliancentury*0.0003032)) %% 360

  ### geom mean anom sun (deg)
  tzpoints_time$geom_anom <- 357.52911 + tzpoints_time$juliancentury*(35999.05029 - 0.0001537*tzpoints_time$juliancentury)

  ### eccent earth orbit
  tzpoints_time$eccent_orbit <- 0.016708634 - tzpoints_time$juliancentury*(0.000042037 + 0.0000001267*tzpoints_time$juliancentury)

  ### sun eq of ctr
  tzpoints_time$sun_eqctr <- sin(tzpoints_time$geom_anom*(pi/180))*(1.914602 - tzpoints_time$juliancentury*(0.004817 + 0.000014*tzpoints_time$juliancentury)) +
    sin((2*tzpoints_time$geom_anom)*(pi/180))*(0.019993 - 0.000101*tzpoints_time$juliancentury) + sin((3*tzpoints_time$geom_anom)*(pi/180))*0.000289

  ### sun true long (deg)
  tzpoints_time$sun_long <- tzpoints_time$geom_long + tzpoints_time$sun_eqctr

  ### sun true anom (deg)
  tzpoints_time$sun_anom <- tzpoints_time$geom_anom + tzpoints_time$sun_eqctr

  ### sun rad vector (AUs)
  tzpoints_time$sun_rad <- (1.000001018*(1 - tzpoints_time$eccent_orbit*tzpoints_time$eccent_orbit))/(1 + tzpoints_time$eccent_orbit*cos(tzpoints_time$sun_anom*(pi/180)))

  ### sun app long (deg)
  tzpoints_time$sun_applong <- tzpoints_time$sun_long - 0.00569-0.00478*sin((125.04 - 1934.136*tzpoints_time$juliancentury)*(pi/180))

  ### mean obliq ecliptic (deg)
  tzpoints_time$obliq_ecliptic <- 23 + (26 + ((21.448 - tzpoints_time$juliancentury*(46.815 + tzpoints_time$juliancentury*(0.00059 - tzpoints_time$juliancentury*0.001813))))/60)/60

  ### obliq corr (deg)
  tzpoints_time$obliq_corr <- tzpoints_time$obliq_ecliptic + 0.00256*cos((125.04 - 1934.136*tzpoints_time$juliancentury)*(pi/180))

  ### sun rt ascen (deg)
  tzpoints_time$sun_rtascen <- (atan2(cos(tzpoints_time$obliq_corr*(pi/180))*sin(tzpoints_time$sun_applong*(pi/180)), cos(tzpoints_time$sun_applong*(pi/180))))*(180/pi)

  ### sun declin (deg)
  tzpoints_time$sun_declin <- (asin(sin(tzpoints_time$obliq_corr*(pi/180))*sin(tzpoints_time$sun_applong*(pi/180))))*(180/pi)

  ### var y
  tzpoints_time$var_y <- tan((tzpoints_time$obliq_corr/2)*(pi/180))*tan((tzpoints_time$obliq_corr/2)*(pi/180))

  ### eq of time (minutes)
  tzpoints_time$eq_time <- 4*((tzpoints_time$var_y*sin(2*(tzpoints_time$geom_long*(pi/180)))-2*tzpoints_time$eccent_orbit*sin(tzpoints_time$geom_anom*(pi/180))+4*tzpoints_time$eccent_orbit*tzpoints_time$var_y*sin(tzpoints_time$geom_anom*(pi/180))*cos(2*(tzpoints_time$geom_long*(pi/180)))-0.5*tzpoints_time$var_y*tzpoints_time$var_y*sin(4*(tzpoints_time$geom_long*(pi/180)))-1.25*tzpoints_time$eccent_orbit*tzpoints_time$eccent_orbit*sin(2*(tzpoints_time$geom_anom*(pi/180))))*(180/pi))

  ### HA sunrise (deg)
  tzpoints_time$ha_sunrise <- (acos(cos(90.833*(pi/180))/(cos(tzpoints_time$Latitude*(pi/180))*cos(tzpoints_time$sun_declin*(pi/180))) - tan(tzpoints_time$Latitude*(pi/180))*tan(tzpoints_time$sun_declin*(pi/180))))*(180/pi)

  ### solar noon (LST)
  tzpoints_time$solar_noon <- (720 - 4 * tzpoints_time$Longitude - tzpoints_time$eq_time + tzpoints_time$timezone*60)/1440
  tzpoints_time$solar_noon_sas <- tzpoints_time$solar_noon*86400 ################################################ double check

  ### accumulate sunrise/sunset times (removing NA values)
  daily_sunrise <- ((tzpoints_time$solar_noon*1440 - tzpoints_time$ha_sunrise*4)/1440)*86400
  daily_sunset  <- ((tzpoints_time$solar_noon*1440 + tzpoints_time$ha_sunrise*4)/1440)*86400
  
  is_valid <- !is.nan(daily_sunrise)
  
  tzpoints_time$sunrise_sum <- ifelse(is_valid,
                                      tzpoints_time$sunrise_sum + daily_sunrise,
                                      tzpoints_time$sunrise_sum)
  tzpoints_time$sunset_sum <- ifelse(is_valid,
                                     tzpoints_time$sunset_sum + daily_sunset,
                                     tzpoints_time$sunset_sum)
  tzpoints_time$valid_days <- tzpoints_time$valid_days + as.integer(is_valid)
  
  ### sunlight duration (minutes)
  tzpoints_time$sunlight_duration <- 8*tzpoints_time$ha_sunrise

  ### true solar time (min)
  tzpoints_time$true_solar_time <- (tzpoints_time$time_clean*1440 + tzpoints_time$eq_time + 4*tzpoints_time$Longitude - 60 * tzpoints_time$timezone) %% 1440

  ### hour angle (deg)
  tzpoints_time$hour_angle <- ifelse((tzpoints_time$true_solar_time / 4) < 0,
                                     tzpoints_time$true_solar_time / 4 + 180,
                                     tzpoints_time$true_solar_time / 4 - 180)

  ### solar zenith angle (deg)
  tzpoints_time$solar_zenith_angle <- (acos(sin(tzpoints_time$Latitude*(pi/180))*sin(tzpoints_time$sun_declin*(pi/180))+cos(tzpoints_time$Latitude*(pi/180))*cos(tzpoints_time$sun_declin*(pi/180))*cos(tzpoints_time$hour_angle*(pi/180))))*(180/pi)

  ### solar elevation angle (deg)
  tzpoints_time$solar_elevation_angle <- 90 - tzpoints_time$solar_zenith_angle

  ### approx atmospheric refraction (deg)
  tzpoints_time$atmospheric_refraction = ifelse(tzpoints_time$solar_elevation_angle > 85, 0,
                                    ifelse((85>=tzpoints_time$solar_elevation_angle) & (tzpoints_time$solar_elevation_angle >5), (58.1/tan(tzpoints_time$solar_elevation_angle*(pi/180)) - 0.07/((tan(tzpoints_time$solar_elevation_angle*(pi/180)))**3) + 0.000086/((tan(tzpoints_time$solar_elevation_angle*(pi/180)))**5))/3600,
                                    ifelse((5 >= tzpoints_time$solar_elevation_angle) & (tzpoints_time$solar_elevation_angle > -0.575), (1735 + tzpoints_time$solar_elevation_angle*(103.4+tzpoints_time$solar_elevation_angle*(-12.79 + tzpoints_time$solar_elevation_angle*0.711)))/3600,
                                    (-20.772 / tan(tzpoints_time$solar_elevation_angle*(pi/180)))/3600)))


  ### solar elevation corrected for atm refraction (deg)
  tzpoints_time$solar_elevation_cor <- tzpoints_time$solar_elevation_angle + tzpoints_time$atmospheric_refraction

  ### solar azimuth angle (deg cv from N)
  tzpoints_time$solar_azimuth <- ifelse(tzpoints_time$hour_angle>0, ((180/pi)*(acos(((sin(tzpoints_time$Latitude*(pi/180))*cos(tzpoints_time$solar_zenith_angle*(pi/180))) - sin(tzpoints_time$sun_declin*(pi/180)))/(cos(tzpoints_time$Latitude*(pi/180))*sin(tzpoints_time$solar_zenith_angle*(pi/180)))))+180) %% 360,
                                        (540-(180/pi)*(acos(((sin(tzpoints_time$Latitude*(pi/180))*cos(tzpoints_time$solar_zenith_angle*(pi/180))) - sin(tzpoints_time$sun_declin*(pi/180)))/(cos(tzpoints_time$Latitude*(pi/180))*sin(tzpoints_time$solar_zenith_angle*(pi/180)))))) %% 360)

}

################################################################################
## Step 3: calculating average of sunrise times and sunset times over 1 year ###
################################################################################

## compute averages from accumulators and keep only what's needed
tzpoints_suntime <- tzpoints_time %>%
  mutate(
    sunrise_avg_tz = sunrise_sum / valid_days,
    sunset_avg_tz  = sunset_sum  / valid_days
  ) %>%
  dplyr::select(Longitude, Latitude, timezone, sunrise_avg_tz, sunset_avg_tz)

cat("\n Sunrise/ sunset time cacluations complete! \n")

## Joining target points with matching border points to calculate suntime differences

#giving target points border values based on time zone
tzpoints_suntime_matches <- tzpoints_suntime %>%
  left_join(
    border_points %>%
      dplyr::select(timezone, sunrise_avg_tz, sunset_avg_tz) %>%
      rename(
        border_sunrise = sunrise_avg_tz,
        border_sunset = sunset_avg_tz
      ),
    by = c("timezone"))

#Calculating sunrise and sunset difference 
valid_matches <- tzpoints_suntime_matches %>%
  # Calculate differences
  mutate(
    sunrise_difference = sunrise_avg_tz - border_sunrise,
    sunset_difference = sunset_avg_tz - border_sunset
  )

#selecting only the variables we need
valid_matches <- valid_matches %>%
  select(Longitude, Latitude, timezone, sunrise_difference, sunset_difference)

# Save final results
cat("\n Saving final datasets \n")

export_complete_path <- paste0(matched_prefix, tile_name, ".rds")

saveRDS(valid_matches, export_complete_path)

cat("\n Tile '", tile_name, "' complete! \n")

