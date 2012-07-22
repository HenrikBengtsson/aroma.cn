dataSet <- "GSE12702"
tags <- "ASCRMAv2"
chipType <- "Mapping250K_Nsp"
ds <- AromaUnitPscnBinarySet$byName(dataSet, tags=tags, chipType=chipType)
print(ds)

# Extract tumors and normals
idxs <- seq(from=1, to=nbrOfFiles(ds), by=2)
dsT <- extract(ds, idxs);
idxs <- seq(from=2, to=nbrOfFiles(ds), by=2)
dsN <- extract(ds, idxs);

# Setup Paired PSCBS model
seg <- PairedPscbsModel(dsT=dsT, dsN=dsN)
print(seg)

# Segment all tumor-normal pairs
fit(seg, verbose=-10)
