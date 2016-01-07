-- -*- coding: utf-8 -*-

newPackage(
        "ResolutionDim2",
        Version => "1.0",
        Date => "January 6, 2016",
        Authors => {{Name => "Guillem Blanco",
                  Email => "gblanco92@gmail.com"}},
        Headline => "...",
        DebuggingMode => false
        )

export {"..."}

load "./puiseuxExpansion.m2";
load "./proximityMatrix.m2";
load "./basePoints.m2";

beginDocumentation()
  load "./docPuiseux.m2";
  load "./docProximity.m2";
  load "./docBase.m2";

TEST
  load "./testPuiseux.m2";
  load "./testProximity.m2";
  load "./testBase.m2";

end
