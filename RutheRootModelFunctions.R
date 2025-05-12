#' Title
#'
#' @param df a data frame with the columns Date, Tsum, TMPM
#' @param rws the relative growth rate parameter of the Richards-Function
#' @param Tbase the base temperature [°C]
#' @param Wsmax the maximum shoot dry matter [g/m²]
#' @param rf the Richards function exponential parameter
#' @param ID a character string for the ID of the treatment
#'
#' @returns a data frame with the columns Date, Tsum, dW_dt, ShootDM
#' @export
#'
# @examples
NumRichards <- function(df, rws, Tbase, Wsmax, rf, ID="") {
  ShootSim <- data.frame(Date=Date(length(length(df$TMPM))), Tsum=numeric(length(df$TMPM)), 
                         dW_dt=numeric(length(length(df$TMPM))),ShootDM=numeric(length(length(df$TMPM))))
  ShootSim[1,"Date"] <- df[1,"Date"]
  ShootSim[1,"Tsum"] <- df[1,"Tsum"]
  ShootSim[1,"ShootDM"] <- 10
  for (i in 2:length(df$TMPM)){
    #  i <- 2
    ShootSim[i,"Date"] <- df[i,"Date"]
    ShootSim[i,"Tsum"] <- df[i,"Tsum"]
    ShootSim[i,"dW_dt"] <- rws*pmax(df[i,"TMPM"]-Tbase,0)*ShootSim[i-1,"ShootDM"] *((Wsmax^rf-ShootSim[i-1,"ShootDM"]^rf)/(rf*Wsmax^rf))
    ShootSim[i,"ShootDM"] <- ShootSim[i-1,"ShootDM"]+ShootSim[i,"dW_dt"] 
  }
  ShootSim$ID <- ID
  return(ShootSim)
}



#' rld_z_t_f
#'
#' @param DMfineRoot # fine root dry matter
#' @param sp_RL # specific root length
#' @param zr # maximum rooting depth  
#' @param Ratio # ratio of RLD on the surface to RLD at the maximum rooting depth
#' @param z1 # upper boundary of the soil layer
#' @param z2 # lower boundary of the soil layer
#'
#' @return average root length density in the soil layer between z1 and z2
#' @export
#'
#' @examples
rld_z_t_f <- function (DMfineRoot=100, sp_RL=7000, zr=120, Ratio=0.01731, z1=40, z2=60)  {

if (z1>zr) {
  WLD_z_t_f <- 0.0
  return(WLD_z_t_f)
  
} else {
  a <- -log(Ratio)/zr
  SRL <- DMfineRoot*sp_RL/1e4
  WLD0 <- (SRL*a)/(1-exp(-a*zr))
  WLD_z_t_f <- WLD0*(exp(-a*z1)-exp(-a*min(z2,zr)))/(a*(z2-z1))
  return(WLD_z_t_f)
}

}



#' RootModDM
#'
#' @param df External data frame with columns for time, temperature, and other parameters
#' @param Tbase Base temperature for root growth
#' @param zrmax Maximum rooting depth
#' @param k_zb Coefficient for root growth rate
#' @param Ratio Ratio of root length density on the surface to root length density at maximum rooting depth
#' @param fFineRoot0 Initial fine root fraction of daily assimilates
#' @param fFineRootDec Decay rate of fine root fraction
#' @param sp_WL Specific root length [cm/g]
#'
#' @returns Data frame with updated root variables
#' @export
#'
# @examples
RootModDM <- function (df, Tbase, zrmax, k_zb, Ratio, fFineRoot0, fFineRootDec, sp_WL ){
  df$zr <- 0
  df$f_root <- 0
  df$dWt_dt <- 0
  df$DM_root <- 0
  df$WL <- 0
  df$DM_root[1] <- 0
  df$rld_15 <- 0
  df$rld_30 <- 0
  df$rld_45 <- 0
  df$rld_60 <- 0
  df$rld_75 <- 0
  df$rld_90 <- 0
  df$rld_105 <- 0
  df$rld_120 <- 0
  for (i in 2:nrow(df)){
    #  i <- 50
    Teff <- max(0,df$TMPM[i]-Tbase)
    df$zr[i] <- min(zrmax,df$zr[i-1] + k_zb*Teff)
    df$f_root[i] <- max(0, fFineRoot0 - fFineRootDec*df$Tsum[i])
    df$dWt_dt[i] <- 1/(1-df$f_root[i])*df$dW_dt[i]
    df$DM_root[i] <- df$DM_root[i-1]+df$dW_dt[i]*df$f_root[i]
    df$WL[i] <- df$DM_root[i]*sp_WL/1e4
    df$rld_15[i]   <- rld_z_t_f(DMfineRoot= df$DM_root[i], sp_WL, df$zr[i], Ratio, 0, 15)
    df$rld_30[i]  <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 15, 30)
    df$rld_45[i]  <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 30, 45)
    df$rld_60[i]  <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 45, 60)
    df$rld_75[i]  <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 60, 75)
    df$rld_90[i]  <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 75, 90)
    df$rld_105[i] <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 90, 105)
    df$rld_120[i] <- rld_z_t_f(DMfineRoot=df$DM_root[i], sp_WL, df$zr[i], Ratio, 105, 120)
  }
  return(df)
}





nls.wrapper.RootModDM <-function(df_ext, Tbase, zrmax, k_zb, Ratio, fFineRoot0, fFineRootDec, sp_WL)
  {

  BVector<- data.frame()
  for (year in unique(df_ext$Year)) {
    for (Nrate in unique(df_ext$Nrate)) {
      # select DM growth  data for the current year and Nrate
#      df_tmp <- df_ext %>% dplyr::filter(Year == year) %>% dplyr::filter( Nrate == Nrate)
      df_tmp <- df_ext[df_ext$Year==year & df_ext$Nrate==Nrate,]
      # model the root growth
      Results <- RootModDM(df_tmp, Tbase, zrmax, k_zb, Ratio, fFineRoot0, fFineRootDec, sp_WL)
      Results <- Results %>%  select(Nrate, Date, rld_15,
                                                                   rld_30,
                                                                   rld_45,
                                                                   rld_60,
                                                                   rld_75,
                                                                   rld_90,
                                                                   rld_105,
                                                                   rld_120)
      
      Results <- Results %>% pivot_longer(cols = c(rld_15,
                                                 rld_30, rld_45, rld_60,
                                                 rld_75, rld_90, rld_105,
                                                 rld_120), names_to = "Depth", values_to = "rld", names_prefix = "rld_")
      
      
      BVector <- rbind(BVector, Results[,c("Date", "Nrate", "Depth","rld")])

    }
  }
  BVector$ID <- paste(BVector$Nrate, BVector$Date, BVector$Depth, sep="_")
  BVector <- as.data.table(BVector)
  setkey(BVector, ID)
#  Results$ObsIDs <- obsIDs
  #BVector <- BVector %>% dplyr::filter (ID %in% obsIDs)
  BVector <- BVector[BVector$ID %in% rld.data$obsIDs,]
  
  #BVector <- BVector %>% arrange(ID)
  return(BVector$rld) # return the root length density values for the specified IDs
}


optim.wrapper.RootModDM <-function(x)
  
{
  # x is a vector of parameters
  Tbase <- 0
  zrmax <- 140
  k_zb <- k_zb
  Ratio <- x[1]
  fFineRoot0 <- x[2]
  fFineRootDec <- x[3]
  sp_WL <- 7000

  BVector<- data.frame()
  for (year in unique(AllShootSim$Year)) {
    for (Nrate in unique(AllShootSim$Nrate)) {
      # select DM growth  data for the current year and Nrate
      #      df_tmp <- AllShootSim %>% dplyr::filter(Year == year) %>% dplyr::filter( Nrate == Nrate)
      df_tmp <- AllShootSim[AllShootSim$Year==year & AllShootSim$Nrate==Nrate,]
      # model the root growth
      Results <- RootModDM(df_tmp, Tbase, zrmax, k_zb, Ratio, fFineRoot0, fFineRootDec, sp_WL)
      Results <- Results %>%  select(Nrate, Date, rld_15,
                                     rld_30,
                                     rld_45,
                                     rld_60,
                                     rld_75,
                                     rld_90,
                                     rld_105,
                                     rld_120)
      
      Results <- Results %>% pivot_longer(cols = c(rld_15,
                                                   rld_30, rld_45, rld_60,
                                                   rld_75, rld_90, rld_105,
                                                   rld_120), names_to = "Depth", values_to = "rld", names_prefix = "rld_")
      
      
      BVector <- rbind(BVector, Results[,c("Date", "Nrate", "Depth","rld")])
      
    }
  }
  BVector$ID <- paste(BVector$Nrate, BVector$Date, BVector$Depth, sep="_")
  BVector <- as.data.table(BVector)
  setkey(BVector, ID)
  #  Results$ObsIDs <- obsIDs
  #BVector <- BVector %>% dplyr::filter (ID %in% obsIDs)
  BVector <- BVector[BVector$ID %in% rld.data$obsIDs,]
  BVector$rld <- ifelse(is.na(BVector$rld) , 0, BVector$rld)
  objective <- sum((log(BVector$rld+1e-4) - log(rld.data$rld))^2) 
  #BVector <- BVector %>% arrange(ID)

  return(objective)
}
