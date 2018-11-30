#' A logistic model function
#'
#' This function allows you to fit a logistic function.
#' @param start_L Maximum length
#' @param x0 initial size
#' @param start_K  carrying capacity
#' @param x0_limit lower limit of the x0 variable or the inflection point of the curve
#' @param near_limit closest range to fit the model or the closest range to a fish observable in both cameras
#' @keywords logisitic model
#' @export
#' @examples
#' logistic_model()
 

logistic_model<-function(species="Fish",range,ob_dens,weights=1,x0_limit=2,near_limit=1){
  require(ggplot2)
  ob_dens<-ob_dens/max(ob_dens)
  data2<-data.frame(range,ob_dens)
  if(sum(data2$ob_dens)>0){
    par1<-c(start_L=-.5,start_K=.1,x0=.01)
    fit1<-optim(par1,logistic_modlike,NULL,range,ob_dens,weights, method="L-BFGS-B",lower=c(-Inf,-Inf,x0_limit))
    fit1
    aic1<-2*3+2*fit1$value
    test1<-data.frame(Group.2=seq(near_limit,6,.1),predicted=logistic_fit(seq(near_limit,6,.1),fit1$par))
    
    density_est<-test1[1,2]
    
    #bootstrap_errors
    boot_error<-array(dim=c(200,1))
    for(i in 1:500){
      data1<-data.frame(range,ob_dens,weights)
      boot1<-data1[sample(1:nrow(data1),replace=TRUE),]
      if(sum(boot1[,2])!=0){
        fit2<-optim(par1,logistic_modlike,NULL,boot1[,1],boot1[,2],boot1[,3], method="L-BFGS-B",lower=c(-Inf,-Inf,x0_limit))
        boot_error[i]<-logistic_fit(1,fit2$par)}
      if(sum(boot1[,2])==0){boot_error[i]<-0}}
    boot_error<-sd(boot_error)
    
    #plotting the output
    p <- ggplot()+geom_point(data=data2, aes(x=range,y=ob_dens)) + ggtitle("Logistic model")+
      geom_line(data=test1, aes(x=Group.2, y=predicted))+geom_text(aes(5, max(ob_dens)*.88, label=paste("AIC =",round(aic1,0))))+xlab("Range bin (m)")+ylab("Volumetric density")+
      geom_text(aes(5, max(ob_dens), label=paste("Start_L =",round(fit1$par[1],4))))+
      geom_text(aes(5, max(ob_dens)*.96, label=paste("Start_k = ",round(fit1$par[2],2))))+
      geom_text(aes(5, max(ob_dens)*.92, label=paste("X0 = ",round(fit1$par[3],2))))
    png(filename=paste(species,"Logistic_Fit.png",sep=""))
    print(p)
    dev.off()
    return(list(start_L=fit1$par[1],start_K=fit1$par[2],x0=fit1$par[3],likelihood=fit1$value,convergence=fit1$convergence,counts=fit1$counts,message=fit1$message,aic=aic1,intcpt=density_est,sdev=boot_error))}
  if(sum(data2$ob_dens)==0){
    return(list(c("no fish",aic=0,intcpt=0,sdev=0)))}
}

#' A logistic fitting function
#'
#' This function fits the logistic function.
#' @param start_L Maximum length
#' @param x0 initial size
#' @param start_K  carrying capacity
#' @keywords logisitic model
#' @export
#' @examples
#' logistic_fit()

logistic_fit<-function(range,par){
  start_L<-par[1]
  start_k<-par[2]
  x0<-par[3]
  out<-start_L/(1+exp(-start_k*(range-x0)))
  return(dens=out)}

#' A negative log-likelihood function
#'
#' This function estimates the likelihood of your parameters
#' @param start_L Maximum length
#' @param x0 initial size
#' @param start_K  carrying capacity
#' @keywords logistic model
#' @export
#' @examples
#' logistic_modlike()

logistic_modlike<-function(par,range,ob_dens,weights=1)
{
  pdens<-logistic_fit(range,par)
  sigma<-sd(ob_dens)
  likj<-ob_dens
  likj<-log(sigma)+0.5*log(2*3.141593)+((ob_dens-pdens)^2)*weights/(2*sigma^2)
  lik_sum=sum(likj)
  return(lik_sum)
}

#' Distance function
#'
#' This function calculates the absolute density from a frame by correcting for
#' the sighting function for the species. 
#' @param start_L output parameter from the logistic model function 
#' @param start_K output parameter from the logistic model function 
#' @param xo output parameter from the logistic model function 
#' @param intercept denisty at x = 0 predicted from the logistic model function 
#' @param range range over which to calculate the density
#' @keywords logistic model
#' @export
#' @examples
#' dist_function()

dist_function<-function(start_L,start_K,x0,intercept,range)
{
  xs<-range
  ys<-logistic_fit(xs,c(start_L,start_K,x0))/intercept
#  print(plot(ys~xs))
  return(ys)
}


#' Absolute density for a frame
#'
#' This function calculates the absolute density from a frame by correcting for
#' the sighting function for the species. 
#' @param start_L output parameter from the logistic model function 
#' @param start_K output parameter from the logistic model function 
#' @param xo output parameter from the logistic model function 
#' @param range_bin series of ranges of densities from a frame;
#' @param density observed densities by range in a frame
#' @keywords logistic model
#' @export
#' @examples
#' frame_density()

frame_density<-function(start_L,start_K,x0,intercept,range_bin,density,near_limit=1){
#  start_L<-species_parameters$start_L[i]
#  start_K<-species_parameters$start_K[i]
#  x0<-species_parameters$x0[i]
#  intercept<-species_parameters$Intercept[i]
#  range_bin<-data5$range_bin
#  density<-data5$density
xs<-seq(near_limit,6,.25)  
ys<-logistic_fit(xs,c(start_L,start_K,x0))/intercept  
limit1<-min(which(ys<.1))
max_range<-xs[limit1]
xys<-data.frame(xs,ys)

d1<-data.frame(range_bin,density)
d1<-merge(d1,xys,by.x="range_bin",by.y="xs",all.x=TRUE)
d2<-d1$density/d1$ys
d2[d1$range_bin>=max_range]<-NA
return(d2)}

#' Zero filled data
#'
#' This function takes the unique deployment, frame, range combinations and 
#' attaches the density estimates for a species where they occur and fills
#' the rest of the unique deployment, frames and range combinations with
#' zeros.. 
#' @param d_data set of density data that includes all unique combinations
#' @param species species of interest
#' @keywords logistic model
#' @export
#' @examples
#' frame_density()
zero_fill<-function(d_data,species){
#  d_data<-density_range
#  species<-"Greenstriped rockfish"
  ddata<-unique(d_data[,1:9])
  sdata<-subset(d_data,d_data$class==species)
  sddata<-merge(ddata,sdata,by=c(colnames(ddata[1:9])),all.x=TRUE)
  sddata$density[is.na(sddata$density)]<-0
  return(sddata)}
