##### Header #####
# Author:       Matthew H. Grinnell
# Affiliation:  Pacific Biological Station, Fisheries and Oceans Canada (DFO)
# Group:        Quantitative Assessment Methods Section, Science
# Address:      3190 Hammond Bay Road, Nanaimo, BC, Canada, V9T 6N7
# Contact:      e-mail: matt.grinnell@dfo-mpo.gc.ca | tel: 250.756.7055
# Project:      Herring
# Code name:    Functions.R
# Version:      1.0
# Date started: Nov 30, 2016
# Date edited:  Jun 02, 2017
#
# Overview:
# Some helper functions that are (or could be) used in more than 1 script.

##### Packages #####
# TODO: There are others!
require(ggplot2)

##### Main #####

# Function to transform from Season to Year
Season2Year <- function(dat) {
  # The herring 'season' column is a combination of the two fishery years: for
  # example, season '20123' indicates the years 2012 and 2013. The input *.dat
  # file for the analysis using ADMB requires this to be an acual year, which we
  # define as the second (i.e., later) year, 2013. This function takes in a
  # vector of seasons (dat), and outputs a vector of years (res).
  # Grab the first 4 characters
  chars <- substr(x = dat, start = 1, stop = 4)
  # Convert to numeric
  digits <- as.numeric(x = chars)
  # Add one to get the year
  res <- digits + 1
  # Return years (as an integer)
  return(as.integer(res))
} # End Season2Year function

# Calculate sum if there are non-NA values, return NA if all values are NA
SumNA <- function(x, omitNA = TRUE) {
  # An alternate version to sum(x, na.rm=TRUE), which returns 0 if x is all NA.
  # This version retuns NA if x is all NA, otherwise it returns the sum.
  # If all NA, NA; otherwise, sum
  ifelse(all(is.na(x)),
    res <- NA,
    res <- sum(x, na.rm = omitNA)
  )
  # Return the result
  return(res)
} # End SumNA function

# Calculate mean if there are non-NA values, return NA if all values are NA
MeanNA <- function(x, omitNA = TRUE) {
  # An alternate version to mean(x, na.rm=TRUE), which returns 0 if x is all NA.
  # This version retuns NA if x is all NA, otherwise it returns the mean.
  # If all NA, NA; otherwise, mean
  ifelse(all(is.na(x)),
    res <- NA,
    res <- mean(x, na.rm = omitNA)
  )
  # Return the result
  return(res)
} # End MeanNA function

# Calculate weighted mean if there are non-NA values, return NA if all values
# are NA
WtMeanNA <- function(x, w, omitNA = TRUE) {
  # An alternate version to weighted.mean(x, w, na.rm=TRUE), which returns 0 if
  # x is all NA. This version retuns NA if x is all NA, otherwise it returns the
  # weighted mean.
  # If all NA, NA; otherwise, weighted mean
  ifelse(all(is.na(x)),
    res <- NA,
    res <- weighted.mean(x, w, na.rm = omitNA)
  )
  # Return the result
  return(res)
} # End MeanNA function

# Calculate maximum if there are non-NA values, return NA if all values are NA
MaxNA <- function(x, omitNA = TRUE) {
  # An alternate version to max(x, na.rm=TRUE), which returns 0 if x is all NA.
  # This version retuns NA if x is all NA, otherwise it returns the maximum.
  # If all NA, NA; otherwise, maximum
  ifelse(all(is.na(x)),
    res <- NA,
    res <- max(x, na.rm = omitNA)
  )
  # Return the result
  return(res)
} # End MaxNA function

# Calculate minimum if there are non-NA values, return NA if all values are NA
MinNA <- function(x, omitNA = TRUE) {
  # An alternate version to min(x, na.rm=TRUE), which returns 0 if x is all NA.
  # This version retuns NA if x is all NA, otherwise it returns the minimum.
  # If all NA, NA; otherwise, minimum
  ifelse(all(is.na(x)),
    res <- NA,
    res <- min(x, na.rm = omitNA)
  )
  # Return the result
  return(res)
} # End MinNA function

# Fill in NA values with rolling mean of previous values
RollMeanNA <- function(dat, n) {
  # Update the NAs in a vector with the mean of the previous values. The number
  # of values used can be less than n for NAs near the start of the vector, and
  # will be a maximum of n for values further along the vector. The value will
  # remain NA if no non-NA values are available.
  # Loop over observations starting with the second observation
  for (i in 2:length(dat)) {
    # If the value is NA
    if (is.na(dat[i])) {
      # Get window for current index: up to n previous values
      muWindow <- (i - min(n, i - 1)):(i - 1)
      # Calculate the mean of the values in the rolling window
      dat[i] <- MeanNA(dat[muWindow])
    } # End if value is NA
  } # End i loop over observations
  # Return the observations with NAs (mostly) replaced by the rolling mean
  return(dat)
} # End RollMeanNA function

# Calculate unique if there are non-NA values, return NA if all values are NA
UniqueNA <- function(x) {
  # An alternate version to unique, which fails sometimes if there are no
  # values. This version retuns NA if x is all NA, otherwise it returns the
  # unique values.
  # If all NA, NA; otherwise, unique
  ifelse(all(is.na(x)),
    res <- NA,
    res <- unique(x)
  )
  # Return the result
  return(res)
} # End UniqueNA function

# Calculate plot margin by expanding the range of x and y as a percentage
CalcXYLims <- function(x, y, pX = 0.01, pY = 0.01) {
  # Set the space between the points and the plot border manually as a specified
  # percentage of the range. Return a list with range for x and y.
  # Get range of x
  rx <- range(x)
  # Calculate difference of x range
  drx <- diff(rx)
  # Calculate new usr for x
  usrX <- c(rx[1] - pX * drx, rx[2] + pX * drx)
  # Get range of y
  ry <- range(y)
  # Calculate difference of y range
  dry <- diff(ry)
  # Calculate new usr for y
  usrY <- c(ry[1] - pY * dry, ry[2] + pY * dry)
  # Return new x and y limits (i.e., xlim and ylim, respectively)
  return(list(x = usrX, y = usrY))
} # End CalcXYLims

# Paste strings nicely
PasteNicely <- function(x, intChars = ", ", nChar = "and ") {
  # Get the length of the vector
  n <- length(x)
  # If there are more than two
  if (n > 2) {
    # Make a print-friendly vector
    x[n] <- paste(nChar, x[n], sep = "")
    # Get print friendly values
    res <- paste(x, collapse = intChars)
  } else { # End if more than two, otherwise
    # Add a space
    nCharSp <- paste(" ", nChar, sep = "")
    # Get print friendly values
    res <- paste(x, collapse = nCharSp)
  } # End if not more than two
  # Return the results
  return(res)
} # End PasteNicely function

# Function to add a new column indicating group ID based on sequential data
ConsecutiveGroup <- function(vec) {
  # Get a vector where the ID depends on whether the value of x is
  # sequential. For example, indicate whether a series of years is sequential,
  # or if there are say three groups of sequential years.
  # Break up the data by groups with consecutive values
  dUniqueGrps <- split(x = vec, f = cumsum(c(1, diff(vec) != 1)))
  # Put group id into the data
  for (g in 1:length(dUniqueGrps)) {
    # Add the group ID to the table
    dUniqueGrps[[g]] <- rep(g, times = length(dUniqueGrps[[g]]))
  } # End g loop over groups
  # Unsplit the list
  GroupID <- as.vector(unlist(dUniqueGrps))
  # Return the data with groups
  return(GroupID)
} # End ConsecutiveGroup function

# Function to add a new column indicating the number of consecutive values
CountConsecutive <- function(vec) {
  # Determine the number of consecutive values in a vector. For example,
  # indicate whether a series of years is sequential, or if there are say three
  # groups of sequential years.
  # Break up the data by groups with consecutive values
  dUniqueGrps <- split(x = vec, f = cumsum(c(1, diff(vec) != 1)))
  # Put group id into the data
  for (g in 1:length(dUniqueGrps)) {
    # Add the group ID to the table
    dUniqueGrps[[g]] <- 1:length(dUniqueGrps[[g]])
  } # End g loop over groups
  # Unsplit the list
  NConsec <- as.vector(unlist(dUniqueGrps))
  # Return the data with groups
  return(NConsec)
} # End ConsecutiveGroup function

# Load shapefiles: land, stat areas, etc
LoadShapefiles <- function(where, a, bMax = 5000) {
  # Load shapefiles for herring sections and aggegate sections to statistical
  # areas and region(s), and make tibbles to plot. In addition, load shapefile
  # for land (i.e., BC coast) and clip to the required extent. Returns a list of
  # shapefiles and tibbles.
  # Message
  cat("Loading shapefiles... ")
  # Get area information
  aSm <- a %>%
    select(SAR, StatArea, Group, Section) %>%
    distinct() %>%
    mutate(
      StatArea = formatC(StatArea, width = 2, flag = "0"),
      Section = formatC(Section, width = 3, flag = "0")
    ) %>%
    arrange(SAR, StatArea, Group, Section)
  # Load the Section shapefile (has Statistical Areas and Regions)
  secRaw <- readOGR(dsn = where$locSec, layer = where$fns$sections, verbose = FALSE)
  # If region is All
  if (region == "All") {
    # Set the region number
    reg_num <- 0
    # Is it a special region?
    reg_type <- "All"
  } else { # End if all, otherwise
    # Get the region number from the table
    reg_num <- regions$SAR[which(regions$Region == region)]
    # Is it a special region?
    reg_type <- regions$Type[which(regions$Region == region)]
  } # End if not all
  # Function to perform some light wrangling
  UpdateSections <- function(dat, keepAll) {
    # Subset the sections to the region(s) in question, and perform some light
    # wrangling to get correct column names and IDs. Note that this would also
    # be the place to apply a slight spatial buffer to ensure that boundaries
    # are contiguous without overlapping.
    # Some light wrangling
    dat@data <- dat@data %>%
      mutate(
        StatArea = as.character(StatArea),
        Section = as.character(Section)
      ) %>%
      select(SAR, StatArea, Section)
    # If retain all the regions
    if (keepAll) {
      # If the region is special
      if (all(reg_type == "Special")) {
        # Update the SAR
        dat@data <- dat@data %>%
          mutate(SAR = ifelse(Section %in% aSm$Section & SAR == -1, 8,
            SAR
          ))
      } # End if special
      # Remove the non-SAR areas
      if (region != "All") res <- dat[dat$SAR != -1, ]
    } else { # End if retain all, otherwise
      # If the region is special
      if (all(reg_type == "Special")) {
        # Subset to the right sections
        res <- dat[dat$Section %in% aSm$Section, ]
        # Update the SAR
        res$SAR <- reg_num
      } else { # End if special, otherwise
        # Subset to the right area
        res <- dat[dat$SAR %in% aSm$SAR, ]
        # Pub sections to subset in proper format
        secSubChar <- formatC(sectionSub, width = 3, format = "d", flag = "0")
        # If requested, get the subset of sections specified
        if (!all(is.na(sectionSub))) res <- res[res$Section %in% secSubChar, ]
      } # End if the region is not special
    } # End if not retaining all
    # Return updated sections
    return(res)
  } # End UpdateSections function
  # Update sections
  secSPDF <- UpdateSections(dat = secRaw, keepAll = FALSE)
  # Convert to data frame and select stat areas in question
  secDF <- secSPDF %>%
    fortify(region = "Section") %>%
    rename(Eastings = long, Northings = lat, Section = group) %>%
    as_tibble()
  # Determine section centroids
  secCent <- gCentroid(spgeom = secSPDF, byid = TRUE)
  # Convert to data frame
  secCentDF <- secCent %>%
    as_tibble() %>%
    rename(Eastings = x, Northings = y) %>%
    mutate(Section = formatC(secSPDF$Section, width = 3, flag = "0")) %>%
    arrange(Section)
  # If 'Groups' has info, dissolve to Groups
  if (!(all(is.na(aSm$Group))) & !region %in% c("JS", "All")) { # & all(is.na(sectionSub))
    # First, remove NAs
    aSmC <- aSm %>%
      filter(!is.na(Group), !is.na(Section)) %>%
      select(StatArea, Group, Section)
    # Merge groups information with sections
    secSPDF@data <- secSPDF@data %>%
      left_join(y = aSmC, by = c("StatArea", "Section"))
    # Dissolve to group
    grpSPDF <- aggregate(x = secSPDF, by = list(Temp = secSPDF$Group), FUN = unique)
    # Convert to data frame
    grpDF <- grpSPDF %>%
      fortify(region = "Group") %>%
      rename(Eastings = long, Northings = lat, Group = group) %>%
      as_tibble()
    # Determine group centroids
    grpCent <- gCentroid(spgeom = grpSPDF, byid = TRUE)
    # Convert to data frame
    grpCentDF <- grpCent %>%
      as_tibble() %>%
      rename(Eastings = x, Northings = y) %>%
      mutate(Group = grpSPDF$Group) %>%
      arrange(Group)
  } else { # End if Groups has info, otherwise
    # No objects (null?)
    grpDF <- NULL
    grpCentDF <- NULL
  } # End if Groups has no info
  # Dissolve to stat area
  saSPDF <- aggregate(x = secSPDF, by = list(Temp = secSPDF$StatArea), FUN = unique)
  # Convert to data frame and select stat areas in question
  saDF <- saSPDF %>%
    fortify(region = "StatArea") %>%
    rename(Eastings = long, Northings = lat, StatArea = group) %>%
    as_tibble()
  # Determine stat area centroids
  saCent <- gCentroid(spgeom = saSPDF, byid = TRUE)
  # Convert to data frame
  saCentDF <- saCent %>%
    as_tibble() %>%
    rename(Eastings = x, Northings = y) %>%
    mutate(StatArea = formatC(saSPDF$StatArea, width = 2, flag = "0")) %>%
    filter(StatArea != "00") %>%
    arrange(StatArea)
  # Dissolve to region
  regSPDF <- aggregate(x = secSPDF, by = list(Temp = secSPDF$SAR), FUN = unique)
  # Convert to data frame and select region(s) in question
  regDF <- regSPDF %>%
    fortify(region = "SAR") %>%
    rename(Eastings = long, Northings = lat, Region = group) %>%
    as_tibble()
  # Get a buffer around the region(s) in question
  buff <- gBuffer(spgeom = regSPDF, width = bMax, byid = FALSE)
  # Calculate the extent
  extBuff <- bbox(buff)
  # Convert the extent to a table
  extDF <- tibble(Eastings = extBuff[1, ], Northings = extBuff[2, ])
  # Determine x:y aspect ratio (for plotting)
  xyRatio <- diff(extDF$Eastings) / diff(extDF$Northings)
  # Read the polygon data: land
  landSPDF <- readOGR(dsn = where$locLand, layer = where$fns$land, verbose = FALSE)
  # Clip the land to the buffer: big
  landCropSPDF <- crop(x = landSPDF, y = extBuff)
  # Convert to data frame
  landCropDF <- landCropSPDF %>%
    fortify(region = "id") %>%
    rename(Eastings = long, Northings = lat) %>%
    as_tibble()
  # If region is All
  if (region == "All") {
    # No need to re-run
    secAllSPDF <- secSPDF
  } else { # End if All, otherwise
    # Update sections (keep all areas)
    secAllSPDF <- UpdateSections(dat = secRaw, keepAll = TRUE)
  } # End if not All
  # Dissolve to stat area
  saAllSPDF <- aggregate(
    x = secAllSPDF, by = list(Temp = secAllSPDF$StatArea),
    FUN = unique
  )
  # Dissolve to region
  regAllSPDF <- aggregate(
    x = secAllSPDF, by = list(Temp = secAllSPDF$SAR),
    FUN = unique
  )
  # Determine region centroids
  regCent <- gCentroid(spgeom = regAllSPDF, byid = TRUE)
  # Get region numbers and names
  temp <- regions %>%
    select(SAR, Region)
  # Convert to data frame
  regCentDF <- regCent %>%
    as_tibble() %>%
    rename(Eastings = x, Northings = y) %>%
    mutate(SAR = regAllSPDF$SAR) %>%
    left_join(y = temp, by = "SAR") %>%
    arrange(SAR)
  # Convert to data frame and select all regions: sections
  secAllDF <- secAllSPDF %>%
    fortify(region = "Section") %>%
    rename(Eastings = long, Northings = lat, Section = group) %>%
    as_tibble()
  # Convert to data frame and select all regions: statistical areas
  saAllDF <- saAllSPDF %>%
    fortify(region = "StatArea") %>%
    rename(Eastings = long, Northings = lat, StatArea = group) %>%
    as_tibble()
  # Convert to data frame and select all regions: regions
  regAllDF <- regAllSPDF %>%
    fortify(region = "SAR") %>%
    rename(Eastings = long, Northings = lat, Region = group) %>%
    as_tibble()
  # Get a buffer around the region(s) in question
  buffAll <- gBuffer(spgeom = regAllSPDF, width = bMax, byid = FALSE)
  # Calculate the extent
  extAllBuff <- bbox(buffAll)
  # Convert the extent to a table
  extAllDF <- tibble(Eastings = extAllBuff[1, ], Northings = extAllBuff[2, ])
  # Determine x:y aspect ration (for plotting)
  xyAllRatio <- diff(extAllDF$Eastings) / diff(extAllDF$Northings)
  # Clip the land to the buffer: big
  landAllCropSPDF <- crop(x = landSPDF, y = extAllBuff)
  # Convert to data frame
  landAllCropDF <- landAllCropSPDF %>%
    fortify(region = "id") %>%
    rename(Eastings = long, Northings = lat) %>%
    as_tibble()
  # Update progress message
  cat("done\n")
  # Return the data frames etc
  return(list(
    secDF = secDF, secCentDF = secCentDF,
    grpDF = grpDF, grpCentDF = grpCentDF,
    saDF = saDF, saCentDF = saCentDF,
    regSPDF = regSPDF, regDF = regDF, regCentDF = regCentDF,
    xyRatio = xyRatio, extDF = extDF,
    landCropSPDF = landCropSPDF, landCropDF = landCropDF,
    secAllSPDF = secAllSPDF,
    secAllDF = secAllDF, saAllDF = saAllDF, regAllDF = regAllDF,
    extAllDF = extAllDF, xyAllRation = xyAllRatio,
    landAllCropDF = landAllCropDF
  ))
} # End LoadShapefiles function

# Function to make a circle
MakeCircle <- function(center = c(0, 0), radius = 1, nPts = 100) {
  # Vector of points
  tt <- seq(from = 0, to = 2 * pi, length.out = nPts)
  # X values (math!)
  xx <- center[1] + radius * cos(tt)
  # Y values (and geometry!)
  yy <- center[2] + radius * sin(tt)
  # Return the data (x and y for a circle)
  return(tibble(X = xx, Y = yy))
} # End MakeCircle function

# Function to switch from 0/1 to No/Yes
YesNo <- function(x) {
  # Input is a column/vector/etc of 0s and 1s, and output is No/Yes as an
  # ordered factor (Yes before No)
  # Update values: 0=No, 1=Yes, otherwise NA
  x <- ifelse(x == 0, "No", ifelse(x == 1, "Yes", x))
  # If if's all Yes/No
  if (all(x %in% c("Yes", "No", NA))) {
    # Make an ordered factor
    xFac <- factor(x, levels = c("Yes", "No"))
  } else { # End if Yes/No, otherwise
    # Just return the data
    xFac <- x
  } # End if not all Yes/No
  # Return updated values
  return(xFac)
} # End YesNo function

# Clip to the extent of a supplied sp object
ClipExtent <- function(dat, spObj, bufDist = NA, silent = FALSE) {
  # Given a set of spatial points, dat, and an spatial polygons object, spObj,
  # determine which points fall inside (or within a given buffer) of the object.
  # Set the spatial X and Y for points that don't overlap (or lie within the
  # buffer) to NA, as they are 'outside' the area (i.e., we assume that the X
  # and Y for these points is wrong) so that we don't show them in charts. We
  # retain the points themselves, because the associated data is still useful.
  # Return the points as a data frame, with updated X and Y info.
  # Require sp library
  require(sp)
  # Creat a buffer around the spatial object, if requested
  if (!is.na(bufDist)) {
    # Require rgeos library
    require(rgeos)
    # Make a buffer
    spObj <- gBuffer(spgeom = spObj, byid = TRUE, width = bufDist)
  } # End if making a buffer
  # Get NAs (if any)
  isNA <- filter(.data = dat, is.na(Eastings), is.na(Northings))
  # Message if there are any
  if (nrow(isNA) > 0 & !silent) {
    cat("Point(s) with missing spatial coordinates (NA):", nrow(isNA), "\n")
  }
  # Wrangle data
  samp <- dat %>%
    filter(!is.na(Eastings) & !is.na(Northings))
  # If there are rows
  if (nrow(samp) > 0) {
    # Make a spatial points object
    spSamp <- SpatialPoints(
      coords = select(samp, Eastings, Northings),
      proj4string = CRS(outCRS)
    )
    # Determine which points are outside the SAR
    inside <- over(x = spSamp, y = spObj)$SAR
    # If any points are outside the SAR
    if (any(is.na(inside))) {
      # Set the X and Y to NA
      samp <- samp %>%
        mutate(
          Eastings = ifelse(is.na(inside), NA, Eastings),
          Northings = ifelse(is.na(inside), NA, Northings)
        )
      # Message
      if (!silent) {
        cat(
          "Point(s) outside SAR boundary: set X and Y to NA:",
          length(inside[is.na(inside)]), "\n"
        )
      }
    } # End if any points are outside the SAR
  } else { # End if there are rows, otherwise
    # Message
    if (!silent) cat("There are no geo-referenced points\n")
  } # End if there are no rows
  # Wrangle data
  res <- samp %>%
    bind_rows(isNA)
  # Return the data
  return(res)
} # End ClipExtent function

# Get the decade (or other rounded value) from the year
GetDecade <- function(dat, r = 10) {
  # Given a vector of years or dates, get the decade, and add an "s" to the end
  # for plot labels. Also deal with NAs.
  # Make sure it's the year
  yr <- as.numeric(format(dat, format = "%Y"))
  # Round to nearest decade (or other value)
  decade <- paste(floor(yr / r) * r, "s", sep = "")
  # Fill in NAs
  decade <- ifelse(decade == "NAs", NA, decade)
  # Return the vector
  return(decade)
} # End GetDecade function

# How to write a long table
WriteLongTable <- function(dat, fn) {
  # Write the xtable (first time)
  print(
    x = dat, file = fn, tabular.environment = "longtable", floating = FALSE,
    include.rownames = FALSE, booktabs = TRUE, only.contents = TRUE,
    NA.string = NA, include.colnames = FALSE, hline.after = FALSE
  )
  # Load the xtable
  xTabLong <- readLines(con = fn, warn = FALSE)
  # Find the midrule
  isMid <- grep(pattern = "midrule", x = xTabLong)
  # Remove the midrule
  xTabLong <- xTabLong[-isMid]
  # Remove the last line if it's empty
  if (xTabLong[length(xTabLong)] == "  ") {
    xTabLong <- xTabLong[-length(xTabLong)]
  }
  # Remove the last line return if there is one
  if (grepl(pattern = "\\ ", x = xTabLong[length(xTabLong)])) {
    xTabLong[length(xTabLong)] <-
      gsub(
        pattern = " \\\\ ", replacement = "", x = xTabLong[length(xTabLong)],
        fixed = TRUE
      )
  }
  # Re-write the xtable
  writeLines(text = xTabLong, con = fn)
  # Add a bottomrule at the end
  #  write( x="\\bottomrule", file=fn, append=TRUE )
} # End WriteLongTable function

# Convert line endings to Linux
ConvertLineEndings <- function(infile) {
  # Grab the text
  txt <- readLines(con = infile)
  # Start a connection (binary)
  f <- file(description = infile, open = "wb")
  # Write file contents
  cat(txt, file = f, sep = "\n")
  # Close the connection
  close(con = f)
} # End ConvertLineEndings function

# Change default ggplot theme to 'black and white'
theme_set(theme_bw())

# Modify default theme
myTheme <- theme(
  legend.box.background = element_rect(fill = alpha("white", 0.7)),
  legend.box.margin = margin(1, 1, 1, 1, "mm"),
  legend.key = element_blank(), legend.margin = margin(), legend.text.align = 1,
  panel.grid.major = element_line(colour = "darkgrey", size = 0.2),
  panel.grid.minor = element_line(colour = "darkgrey", size = 0.1),
  legend.background = element_rect(fill = "transparent"),
  # panel.spacing.x=unit(3, "lines"),
  plot.margin = unit(c(0.1, 0.6, 0.1, 0.1), "lines")
)

# Change numbers into words (1:9)
Num2Word <- function(x) {
  # Get the list of numbers and words
  vec <- c(
    one = 1, two = 2, three = 3, four = 4, five = 5, six = 6, seven = 7, eight = 8,
    nine = 9
  )
  # Get the index corresponding to the number
  ind <- which(vec == x)
  # Get the name
  res <- names(vec[ind])
  # Error if the name isn't there
  if (length(res) == 0) stop("Numbers can be from 1 to 9 only.")
  # Return the result
  return(res)
} # End function Num2Word

# Install missing packages and load required packages (if required)
UsePackages <- function(pkgs, locn = "https://cran.rstudio.com/") {
  # Reverse the list
  rPkgs <- rev(pkgs)
  # Identify missing (i.e., not yet installed) packages
  newPkgs <- rPkgs[!(rPkgs %in% installed.packages()[, "Package"])]
  # Install missing packages if required
  if (length(newPkgs)) install.packages(newPkgs, repos = locn)
  # Loop over all packages
  for (i in 1:length(rPkgs)) {
    # Load required packages using 'library'
    eval(parse(text = paste("suppressPackageStartupMessages(library(", rPkgs[i],
      "))",
      sep = ""
    )))
  } # End i loop over package names
} # End UsePackages function

# Calculate percent change, difference, etc
DeltaPercent <- function(x, nYrs = 1, type) {
  # Numerator: difference
  top <- x - lag(x, n = nYrs)
  # Denominator: depends on the type
  if (type == "PctChange") bot <- lag(x, n = nYrs) # Previous value
  if (type == "PctDiff") bot <- (x + lag(x, n = nYrs)) / 2 # Mean of values
  # Error if there's no denominator
  if (!exists("bot")) stop("No denominator; check 'type'")
  # Calculate percent change
  res <- top / bot * 100
  # Return the result
  return(res)
} # End DeltaPercent function

# Make sure area info is consistent using a spatial overlay
CheckSpatialOverlay <- function(pts, shape, type) {
  # If there are NAs in spatial info
  if (any(is.na(pts$Eastings), is.na(pts$Northings))) {
    # Message re NAs
    warning("NAs in spatial info: Eastings (", sum(is.na(pts$Eastings)),
      "), and Northings (", sum(is.na(pts$Northings)), ")",
      sep = ""
    )
    # Remove the NAs
    pts <- pts %>%
      filter(!is.na(Eastings), !is.na(Northings))
  } # End if there are NAs
  # Subset the points and get stat area and section from the Locations table
  pts <- pts %>%
    rename(StatAreaLoc = StatArea, SectionLoc = Section) %>%
    mutate(Source = type) %>%
    select(
      Source, StatAreaLoc, SectionLoc, LocationCode, LocationName,
      Eastings, Northings, Longitude, Latitude
    ) %>%
    distinct()
  # Convert to spatial object
  ptsSPDF <- SpatialPointsDataFrame(
    coords = select(pts, Eastings, Northings), proj4string = CRS(outCRS),
    data = pts
  )
  # Spatial overlay
  overPts <- over(x = ptsSPDF, y = shape) %>%
    as_tibble() %>%
    mutate(StatArea = as.integer(StatArea), Section = as.integer(Section)) %>%
    rename(StatAreaPoly = StatArea, SectionPoly = Section) %>%
    select(StatAreaPoly, SectionPoly)
  # Bind the spatial info
  res <- bind_cols(pts, overPts) %>%
    filter(StatAreaLoc != StatAreaPoly | SectionLoc != SectionPoly)
  # Warning if there are oddballs
  if (nrow(res) >= 1) {
    warning("Spatial overlay inconsistencies for ", nrow(res), " points")
  }
  # Return the result
  return(res)
} # End CheckAreas function

##### End #####

# Print end of file message
cat("Loaded helper functions: 'Functions.R'\n")
