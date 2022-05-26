---
layout: default
title: ablate::domain::modifiers::Modifier
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::domain::modifiers::Modifier](./ablate::domain::modifiers::Modifier.html)
## ablate::domain::modifiers::DistributeWithGhostCells
Distribute DMPlex with ghost cells

ghostCellDepth
: (int) the number of ghost cells to share on the boundary.  Default is 1.

## ablate::domain::modifiers::GhostBoundaryCells
Adds ghost cells to the boundary

labelName
: (string) The label specifying the boundary faces, or "Face Sets" if not specified

## ablate::domain::modifiers::SetFromOptions
Sets the specified options on the dm.

## ablate::domain::modifiers::CreateLabel
Creates a new label for all positive points in the function

region (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region describing the new label

function (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the function to evaluate

depth
: (int) The depth in which to apply the label.  The default is zero or cell/element

## ablate::domain::modifiers::TagLabelBoundary
Creates a new label for all faces on the outside of the boundary

region (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to tag the boundary

boundaryFaceRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the new region for the boundary faces

boundaryCellRegion
: ([ablate::domain::Region](./ablate::domain::Region.html)) the new region for the boundary cells

## ablate::domain::modifiers::MergeLabels
Creates a new label for all faces on the outside of the boundary

mergedRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the merged region to create

regions (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html) list) the regions to include in the new merged region

## ablate::domain::modifiers::IntersectLabels
Creates a new label that intersections the provided regions

intersectRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the intersect region to create/used

regions (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html) list) the regions to include in the intersection

## ablate::domain::modifiers::SubtractLabel
Cuts/removes the given region (difference = minuend - subtrahend)

differenceRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the result of the operation

minuendRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the minuend region

subtrahendRegions (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html) list) the region(s) to be removed

incompleteLabel
: (bool) determines if the DMPlexLabelComplete function for the new label is called. (true = DMPlexLabelComplete not called, false = DMPlexLabelComplete called, default is false

## ablate::domain::modifiers::TagLabelInterface
Class to label/tag all faces/cells on the interface between two labels.  The left/right designations are just used to separate the left/right labels.

leftRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the "left" region

rightRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the "right" region

boundaryFaceRegion (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the new region for the newly tagged boundary faces

leftBoundaryCellRegion
: ([ablate::domain::Region](./ablate::domain::Region.html)) optional new region to tag the boundary cells on the "left" of region

rightBoundaryCellRegion
: ([ablate::domain::Region](./ablate::domain::Region.html)) optional new region to tag the boundary cells on the "right" of region

