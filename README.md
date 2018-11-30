Density Calculation Functions for Stationary Platform
================
Chris Rooper
November 29, 2018

Purpose
-------

The functions contained in this package are used to estimate densities of fish using a stationary camera platform (e.g. Williams et al. 2014). The functions follow the method of Williams et al. (2018) In short, this analysis takes a matrix of densities from individual camera frames for a species in range bins with volume estimates for each range bin. It makes a volume-weighted estimate of the detection function by fitting a logistic function. This function can then be used to estimate a density for the species in each frame. These can be averaged across frames to come up with a volumetric density estimate.

![](TrigCamDensityEstimation_files/figure-markdown_github/Figure%201-1.png)

Package Functions and Example Using Bocaccio Rockfish Data
----------------------------------------------------------

As an example, we will use bocaccio rockfish data collected at 30 sites over ~24 hour deployments on the Footprint Bank. These data were collected in October 2017 aboard the R/V Bell-Shimada. In this case the columns are a deployment number, a range bin and a volumetric density in that range bin.

``` r
data(bocaccio)
pander::pandoc.table(bocaccio[1:6,], digits=3, split.tables=Inf)
```

    ## 
    ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    ##    deployment_ID     frame_number   range_bin   unit   site         frame_time          primary_habitat_type   secondary_habitat_type   volume   class   density 
    ## ------------------- -------------- ----------- ------ ------ ------------------------- ---------------------- ------------------------ -------- ------- ---------
    ##  D20171013-T024658        24          1.25       2      4     2017-10-13 02:58:31.276       Low Bedrock               Boulder           0.461     NA        0    
    ## 
    ##  D20171013-T024658        24          1.75       2      4     2017-10-13 02:58:31.276       Low Bedrock               Boulder           0.892     NA        0    
    ## 
    ##  D20171013-T024658        24          2.25       2      4     2017-10-13 02:58:31.276       Low Bedrock               Boulder            1.27     NA        0    
    ## 
    ##  D20171013-T024658        24          2.75       2      4     2017-10-13 02:58:31.276       Low Bedrock               Boulder            1.71     NA        0    
    ## 
    ##  D20171013-T024658        24          3.25       2      4     2017-10-13 02:58:31.276       Low Bedrock               Boulder            2.02     NA        0    
    ## 
    ##  D20171013-T024658        24          3.75       2      4     2017-10-13 02:58:31.276       Low Bedrock               Boulder            2.43     NA        0    
    ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------

To fit the model, we call the logistic model function that serves as a wrapper for the model fitting. The logistic model function is the wrapper that inputs the observations of density, volume and range and fits the logistic model. It also takes inputs of a species name and a deployment number (for ease of looping across species and deployments). Weights can also be included for weighting the AIC calculation. It outputs an estimate the y-intercept of the model (which is used to scale the results)) and a plot of the fitted logistic function (output in ggplot). The output also includes a list of fitted values (parameters, AIC, etc.).

``` r
f1<-logistic_model("Bocaccio rockfish",bocaccio$range_bin,bocaccio$density,bocaccio$volume,2,.5)
print(f1)
```

    ## $start_L
    ##  start_L 
    ## 1.382278 
    ## 
    ## $start_K
    ##   start_K 
    ## -1.278837 
    ## 
    ## $x0
    ## x0 
    ##  2 
    ## 
    ## $likelihood
    ## [1] -1.047705
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $counts
    ## function gradient 
    ##       11       11 
    ## 
    ## $message
    ## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
    ## 
    ## $aic
    ## [1] 3.904589
    ## 
    ## $intcpt
    ## [1] 1.205269

This function will produce a figure of the detecability function that is scaled to 1 at the nearest volume bin (in this case 0.5 m)

![](TrigCamDensityEstimation_files/figure-markdown_github/Figure%202-1.png)

Here is the code for the logistic model function. It calls two other functions; logistic\_fit and logistic\_modlike which make the predictions and estimate the negative log-likelihood respectively.

``` r
logistic_model
```

    ## function (species = "Fish", range, ob_dens, weights = 1, x0_limit = 2, 
    ##     near_limit = 1) 
    ## {
    ##     ob_dens <- aggregate(ob_dens, by = list(range), FUN = "sum")
    ##     if (length(weights == 1)) {
    ##         weights <- rep(1, length(range))
    ##     }
    ##     weights <- aggregate(weights, by = list(range), FUN = "mean")
    ##     ob_dens$x <- ob_dens$x/max(ob_dens$x)
    ##     data2 <- data.frame(range = ob_dens$Group.1, ob_dens = ob_dens$x)
    ##     if (sum(data2$ob_dens) > 0) {
    ##         par1 <- c(start_L = -0.5, start_K = 0.1, x0 = 0.01)
    ##         fit1 <- optim(par1, logistic_modlike, NULL, data2$range, 
    ##             data2$ob_dens, weights$x, method = "L-BFGS-B", lower = c(-Inf, 
    ##                 -Inf, x0_limit))
    ##         fit1
    ##         aic1 <- 2 * 3 + 2 * fit1$value
    ##         test1 <- data.frame(Group.2 = seq(near_limit, 6, 0.1), 
    ##             predicted = logistic_fit(seq(near_limit, 6, 0.1), 
    ##                 fit1$par))
    ##         density_est <- test1[1, 2]
    ##         p <- ggplot() + geom_point(data = data2, aes(x = range, 
    ##             y = ob_dens)) + ggtitle("Logistic model") + geom_line(data = test1, 
    ##             aes(x = Group.2, y = predicted)) + geom_text(aes(5, 
    ##             max(data2$ob_dens) * 0.88, label = paste("AIC =", 
    ##                 round(aic1, 0)))) + xlab("Range bin (m)") + ylab("Volumetric density") + 
    ##             geom_text(aes(5, max(data2$ob_dens), label = paste("Start_L =", 
    ##                 round(fit1$par[1], 4)))) + geom_text(aes(5, max(data2$ob_dens) * 
    ##             0.96, label = paste("Start_k = ", round(fit1$par[2], 
    ##             2)))) + geom_text(aes(5, max(data2$ob_dens) * 0.92, 
    ##             label = paste("X0 = ", round(fit1$par[3], 2))))
    ##         png(filename = paste(species, "_Logistic_Fit.png", sep = ""))
    ##         print(p)
    ##         dev.off()
    ##         return(list(start_L = fit1$par[1], start_K = fit1$par[2], 
    ##             x0 = fit1$par[3], likelihood = fit1$value, convergence = fit1$convergence, 
    ##             counts = fit1$counts, message = fit1$message, aic = aic1, 
    ##             intcpt = density_est))
    ##     }
    ##     if (sum(data2$ob_dens) == 0) {
    ##         return(list(c("no fish", aic = 0, intcpt = 0)))
    ##     }
    ## }
    ## <bytecode: 0x000000000bade238>
    ## <environment: namespace:TrigCamDensityEstimation>

The logistic\_fit function takes the parameter values and calculates a predicted density based on the range. This is a logistic function (e.g. logistic growht model) The parameters from this equation are fit using negative-log-likelihood estimation (see logistic\_modlike function below).

``` r
logistic_fit
```

    ## function (range, par) 
    ## {
    ##     start_L <- par[1]
    ##     start_k <- par[2]
    ##     x0 <- par[3]
    ##     out <- start_L/(1 + exp(-start_k * (range - x0)))
    ##     return(dens = out)
    ## }
    ## <bytecode: 0x000000000cf1a498>
    ## <environment: namespace:TrigCamDensityEstimation>

The logistic\_modlike function computes the negative-log-likelihood for the model based on parameter estimates (this is the function you are minimizing). For this function you can include optional weights, if for example you want to down-weight the points closest to the camera with low volume, but potentially high density. Default weights are = 1

``` r
logistic_modlike
```

    ## function (par, range, ob_dens, weights = 1) 
    ## {
    ##     pdens <- logistic_fit(range, par)
    ##     sigma <- sd(ob_dens)
    ##     likj <- ob_dens
    ##     likj <- log(sigma) + 0.5 * log(2 * 3.141593) + ((ob_dens - 
    ##         pdens)^2) * weights/(2 * sigma^2)
    ##     lik_sum = sum(likj)
    ##     return(lik_sum)
    ## }
    ## <bytecode: 0x000000000cf1c8f8>
    ## <environment: namespace:TrigCamDensityEstimation>

Below is an example implementation which can be looped for multiple species and combining into a Footprint-wide abundance estimate. This example is how I calculated the biomass estimates for the UHSI project (Rooper et al. in prep). Basically, it takes the detecability function from above and corrects the densities by frame for the distance range, then averages the densities and computes out an overall density. This is then expanded to an areal estimate using the Footprint Bank area estimate for rocky habitat.

Example to Generate Abundance Estimate
--------------------------------------

``` r
#######DATA ANALYSES####################################
species_list<-c("Bocaccio rockfish")
deployment<-unique(bocaccio$deployment_ID)
species_parameters<-array(dim=c(0,5))
colnames(species_parameters)<-c("Species","start_L","start_K","x0","Intercept")


#fit the logistic model to density by range data for each species
for(i in 1:length(species_list)){
f1<-logistic_model(species_list[i],bocaccio$range_bin,bocaccio$density,bocaccio$volume,2,.5)
temp1<-data.frame(Species=species_list[i],start_L=f1[1],start_K=f1[2],x0=f1[3],Intercept=f1[9])
species_parameters<-rbind(species_parameters,temp1)
}

colnames(species_parameters)<-c("Species","start_L","start_K","x0","Intercept")
functs<-array(dim=0,3)
Range<-seq(.5,6,.1)
for(i in 1:length(species_list)){
x1<-dist_function(species_parameters$start_L[i],species_parameters$start_K[i],species_parameters$x0[i],species_parameters$Intercept[i],Range)
x1<-data.frame(Range,x1,species_list[i])
functs<-rbind(functs,x1)
}
colnames(functs)<-c("Range","Relative_density","Species")

dev.new(width=9,height=3,units="in",noRStudioGD=TRUE)
p<-ggplot(functs)+geom_line(aes(y=Relative_density,x=Range,colour=Species))+xlab("Range (m)")+
  ylab("Relative density")+geom_hline(yintercept=0.10,linetype="dashed",color="red",size=.75)
print(p)

density_by_frame<-array(dim=c(0,4))
colnames(density_by_frame)<-c("Group.1","Group.2","x","species")
for(i in 1:length(species_list)){
#  data5<-zero_fill(density_range,species_list[i])
  density_by_frame1<-aggregate(bocaccio$density,by=list(bocaccio$deployment_ID,bocaccio$frame_number,bocaccio$range_bin),FUN="sum")
  density_by_frame1$corrected_density<-frame_density(species_parameters$start_L[i],species_parameters$start_K[i],species_parameters$x0[i],species_parameters$Intercept[i],density_by_frame1$Group.3,density_by_frame1$x)
  density_by_frame1<-aggregate(density_by_frame1$corrected_density,by=list(density_by_frame1$Group.1,density_by_frame1$Group.2),FUN="mean",na.rm=TRUE)
  density_by_frame1$species<-species_list[i]
  density_by_frame<-rbind(density_by_frame,density_by_frame1)
}  
colnames(density_by_frame)<-c("deployment_ID","frame_number","density","species")  

#ABUNDANCE ESTIMATION
area<-10*10*6819*1.25

#USING ALL TIME
fd1<-aggregate(density_by_frame$density,by=list(density_by_frame$species,density_by_frame$deployment_ID),FUN="mean")
fd2<-aggregate(fd1$x,by=list(fd1$Group.1),FUN="mean")
fd3<-aggregate(fd1$x,by=list(fd1$Group.1),FUN="sd")
abundance_estimate<-merge(fd2,fd3,by="Group.1")
colnames(abundance_estimate)<-c("Species","Density","Density_SD")
abundance_estimate$Population<-abundance_estimate$Density*area
abundance_estimate$Population_SE<-sqrt(abundance_estimate$Density_SD^2*(area^2))/sqrt(30)
abundance_estimate$CV<-abundance_estimate$Density_SD/abundance_estimate$Density
pandoc.table(abundance_estimate)
```

    ## 
    ## --------------------------------------------------------------------------------
    ##       Species        Density    Density_SD   Population   Population_SE    CV   
    ## ------------------- ---------- ------------ ------------ --------------- -------
    ##  Bocaccio rockfish   0.001307    0.002245       1114          349.4       1.718 
    ## --------------------------------------------------------------------------------

References
----------

Rooper CN, Williams K, Towler R. in prep. Estimating rockfish habitat-specific abundance and behaviour using stationary cameras in the southern California Bight.

Williams K, De Robertis A, Berkowitz Z, Rooper C, Towler R. 2014. An underwater stereo-camera trap. Methods in Oceanography 11:1-12. <http://dx.doi.org/10.1016/j.mio.2015.01.003>

Williams K, Rooper CN, De Robertis A, Levine M, Towler R. 2018. A method for computing volumetric fish density using stereo cameras. Journal of Experimental Marine Biology and Ecology 508:21-26. <https://doi.org/10.1016/j.jembe.2018.08.001>
