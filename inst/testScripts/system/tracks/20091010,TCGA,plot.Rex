if (interactive()) savehistory();
library("aroma.affymetrix");
library("aroma.cn");
library("R.menu");

# Find the directory of this script being source():ed.
path <- getParent(names(findSourceTraceback())[1]);
pathname <- Arguments$getReadablePathname("utils.R", path=path);
source(pathname);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figRootPath <- Arguments$getWritablePath("figures/talk/");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

loadData <- function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  rootPath <- Arguments$getReadablePath("totalAndFracBData");
  dataSets <- list.files(path=rootPath);
  dataSet <- textMenu(dataSets, value=TRUE);
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load total DNA signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dsT <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*");
  print(dsT);
  
  # Assert that annotation data exists
  ugp <- getAromaUgpFile(dsT);
  print(ugp);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load allele B fractions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dsB <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*");
  print(dsB);
  
  dsList <- list(total=dsT, fracB=dsB);
  
  if (getName(dsList$total) == "TCGA") {
    dsList <- lapply(dsList, FUN=extractTcgaPairs);
  }
  
  print(dsList);
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select sample
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ds <- dsList$total;
  sampleNames <- unique(getNames(ds));
  sampleName <- textMenu(sampleNames, value=TRUE);
  
  dsList <- lapply(dsList, FUN=function(ds) {
    extract(ds, sampleName);
  });
  print(dsList);
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  chromosomes <- getChromosomes(ugp);
  chromosome <- textMenu(chromosomes, value=TRUE);
  chrTag <- sprintf("Chr%02d", chromosome);
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- extractTcgaSignals(dsList, chromosome=chromosome, ...);

  attr(data, "dsList") <- dsList;

  print(data);
  data;
} # loadData();


main <- function(data, pa=NULL, region=NULL) {
  if (is.null(pa)) {
    pa <- textPlotAnnotations(pa="default");
  }

  dsList <- attr(data, "dsList");
  ds <- dsList[[1]];
  ugp <- getAromaUgpFile(ds);
  platform <- getPlatform(ds);
  chipType <- getChipType(ds);
  technologyTag <- paste(c(platform, chipType), collapse=",");
  print(technologyTag);
  figPath <- filePath(figRootPath, technologyTag);
  figPath <- Arguments$getWritablePath(figPath);
  print(figPath);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select tracks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  tracksDefault <- names(data);
  tracksDefault[is.element(tracksDefault, c("muN", "rhoN", "rhoT", "rhoTN"))] <- NA;
  tracksDefault <- tracksDefault[!is.na(tracksDefault)];
  tracks <- tracksDefault;
  #selected <- is.element(names(data), tracksDefault);
  #tracks <- selectMenu(names(data), selected=selected);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  attachLocally(pa);
  
  chromosome <- data$C$chromosome;
  dataT <- NULL;

  plot <- "";
  
  while (TRUE) {
    if (is.null(dataT)) {
      dataT <- lapply(data[tracks], FUN=function(obj) {
        if (!is.null(region) && inherits(obj, "RawGenomicSignals"))
          obj <- extractRegion(obj, region=region);
        obj;
      });
      print(dataT);
    }
  
    if (smoothCol != "NA") {
      smoothTracks <- grep("^(C|rho.*Het)$", names(dataT), value=TRUE);
      for (key in smoothTracks) {
        gs <- dataT[[key]];
        if (is.null(gs$smooth)) {
          gs$smooth <- binnedSmoothing(gs, by=200e3);
          dataT[[key]]$gs;
        }
      }
    }
  
    if (plot != "") {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Plot tracks
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      sampleName <- data$C$name;
      chrTag <- sprintf("Chr%02d", as.integer(chromosome));
      if (is.null(region)) {
        regionTag <- NULL;
      } else {
        regionTag <- sprintf("%.2f-%.2fMb", region[1]/1e6, region[2]/1e6);
      }
      paTag <- digest(pa, algo="crc32");

      for (kk in seq(along=dataT)) {
        key <- names(dataT)[kk];
        trackTag <- key;
        tags <- c(chrTag, regionTag, trackTag);
        figName <- paste(c(sampleName, tags), collapse=",");
        print(figName);
  
        gs <- dataT[[kk]];
        if (plot == "png") {
          filename <- sprintf("%s,%s.png", figName, paTag);
          pathname <- file.path(figPath, filename);
          devNew("png", pathname, width=width, height=height);
        } else if (plot == "sceen") {
          devSet(figName);
        }
  
        par(mar=c(3.5,3.5,1.3,1.1), mgp=c(1.8, 0.8, 0));
        par(cex.axis=1.5, cex.lab=1.5);
        par(las=1); # Horizontal numbering of axes
        ylim <- gs$ylim;
        yAt <- gs$yAt;
        axes <- (is.null(yAt));
        if (smoothCol != "NA") {
          col <- rawCol;
          plot(gs, col=col, xlab=xlab, ylab=ylab, ylim=ylim, xScale=xScale, axes=axes);
          gsS <- gs$smooth;
          points(gsS, cex=1.5*smoothCex, col="white", xScale=xScale);
          points(gsS, cex=smoothCex, col=smoothCol, xScale=xScale);
        } else {
          col <- gs$col;
          if (is.null(col)) {
            col <- rawCol;
          }
          plot(gs, col=col, xlab=xlab, ylab=ylab, ylim=ylim, xScale=xScale, axes=axes);
        }
        if (!axes) {
          box();
          axis(side=1);
          axis(side=2, at=yAt);
        }
        if (showSampleName) {
          stext(side=3, pos=0, sampleName, cex=1.5);
        }
        if (showChr) {
          stext(side=3, pos=1, chrTag, cex=1.5);
        }
  #      stext(side=4, pos=0, trackTag, cex=0.8, col="#666666");
  
        devDone();
      } # for (kk ...)
      plot <- "";
    }
  
    choices <- c("Change chromosome", "Change region", "Change tracks", "Settings", "Plot on screen", "Write PNG", "Close all figures");
    ans <- textMenu(choices);
    if (ans == 1) {
      chromosome0 <- chromosome;
      chromosomes <- getChromosomes(ugp);
      chromosome <- textMenu(chromosomes, value=TRUE);
      chrTag <- sprintf("Chr%02d", chromosome);
      data <- extractTcgaSignals(dsList, chromosome=chromosome);
      dataT <- NULL;
      region <- NULL;
    } else if (ans == 2) {
      if (!is.null(region)) {
        region <- region / 1e6;
      }
      region0 <- region;
      region <- textRegion(region=region);
      region <- 1e6 * region;
      if (!isTRUE(all.equal(region, region0))) {
        dataT <- NULL;
      }
    } else if (ans == 3) {
      selected <- is.element(names(data), tracks);
      tracks0 <- tracks;
      tracks <- selectMenu(names(data), selected=selected);
      if (!isTRUE(all.equal(tracks, tracks0))) {
        dataT <- NULL;
      }
    } else if (ans == 4) {
      pa0 <- pa;
      pa <- textPlotAnnotations(pa);
      if (!isTRUE(all.equal(pa, pa0))) {
        attachLocally(pa);
        str(pa);
        graphics.off();
      }
    } else if (ans == 5) {
      plot <- "screen";
    } else if (ans == 6) {
      plot <- "png";
    } else if (ans == 7) {
      graphics.off();
    }
  } # while()
} # main()

if (!exists("pa", mode="list")) {
  pa <- textPlotAnnotations(pa="default");
}

if (!exists("data", mode="list")) {
  data <- loadData(kappa=0.68);
  data <- loadData(kappa=0.58);
}

main(data=data, pa=pa);
