-- -*- coding: utf-8 -*-

newPackage(
        "ResolutionDim2",
        Version => "1.0",
        Date => "January 13, 2016",
        Authors => {{Name => "Guillem Blanco",
                  Email => "gblanco92@gmail.com"}},
        Headline => "Computes the log-resolution of any two-dimensional ideal.",
        DebuggingMode => false,
        AuxiliaryFiles => true
        )

export {"PuiseuxSeries","puiseuxExpansion","Bits","Terms","proximityMatrix","ExtraPoint","basePoints"}

load "./ResolutionDim2/puiseuxSeries.m2";
load "./ResolutionDim2/puiseuxExpansion.m2";
load "./ResolutionDim2/proximityMatrix.m2";
load "./ResolutionDim2/basePoints.m2";

beginDocumentation()
  load "./ResolutionDim2/docPuiseux.m2";
  load "./ResolutionDim2/docProximity.m2";
  load "./ResolutionDim2/docBase.m2";

TEST
  load "./ResolutionDim2/testPuiseux.m2";
  load "./ResolutionDim2/testProximity.m2";
  load "./ResolutionDim2/testBase.m2";

end
