Empirical data used in the manuscript titled: "When traits vary across species but performance doesn't: One solution to the paradox".

The full dataset was first presented in Nolting et al. (2021). Full citation:

Kristen M Nolting, Rachel Prunier, Guy F Midgley, Kent E Holsinger, Intraspecific trait variation influences physiological performance and fitness in the South Africa shrub genus Protea (Proteaceae), Annals of Botany, Volume 127, Issue 4, 1 April 2021, Pages 519–531, https://doi.org/10.1093/aob/mcaa060

For the present analyses we subset the data to 132 observations, and included eight structural traits and six physiological performance traits. We also included site and species IDs as identifying information. Below are descriptions of the column names in the Protea_data.R dataframe (after running through the Data_Wrangling.R script in the \Code folder), which is used in all analyses (excluding the simulation).

1. Species: species identification (PRPU = Protea punctata; PRRE = P. repens; PRLA = P. laurifolia; PRNI = P. nitida; PREX = P. eximia)
2. Site: site sampled

All following are traits. Note that "_scaled" indicates the raw values were centered to a mean of zero and scaled to a standard deviation of one, with respect to all observations. See Nolting et al. 2021 for more detailed information on trait sampling and measurements.

3. BW_scaled: bark width/thickness (mm)
4. WD_scaled: wood density (g cm^-3)
5: Leaf_Area_scaled: area of a single side of the leaf (cm^2)
6: LMA_scaled: leaf mass per area (g cm^-2)
7: LD_scaled: lamina density (g cm^-3)
8: LWR_scaled:  lamina length-width ratio (unitless)
9: Stom_L_scaled: stomatal length (mm)
10: Stom_D_scaled: stomatal density on top leaf surface (mm^-2)
11: Ks_scaled: stem-specific hydraulic conductance
12: LSC_scaled: leaf-specific conductivity
13: Photo_scaled: area-based light-saturated photosynthetic rate 
14: Total_Assim_scaled: leaf-specific photosynthetic rate
15: Cond_scaled: stomatal conductance
16: WUE_Instan_scaled: instantaneous water-use efficiency

