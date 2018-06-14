require(ggplot2)
require(SDMTools)
require(MASS)
require(dplyr)

#Download the csv file here -> 
#Once downloaded, run script & compare fwds or dmen using the FVF and DVD functions
#Examples... FVF("SIDNEY.CROSBY","CONNOR.MCDAVID")
#             DVD("ERIK.KARLSSON","MIKHAIL.SERGACHEV")

  #Load data
    
    Data = read.csv(file.choose(), header = T)
  
  #Filter latest 3-year span
    
    Data_L3 = Data %>% 
                filter(Latest.Season == 2018)
    
  #Filter by position and compute raw scoring rates
    
    Fwds = Data %>%
            filter(pos == "F" & GP >= 41) %>%
            mutate(Rate = PP / TOI * 60)
    
    Dmen = Data %>% 
            filter(pos == "D" & GP >= 41) %>%
            mutate(Rate = PP / TOI * 60)
    
    #Adjust scoring rates and replace zeros with fudge factor
    
    #Compute Conversion Rates
    
      conversion.rates.F = Fwds %>%
                          group_by(Latest.Season) %>%
                          summarise(Mean = mean(Rate))
      
      mean_L3.F = conversion.rates.F %>%
                  filter(Latest.Season == 2018) %>%
                  .$Mean
      
      Fwds = conversion.rates.F %>%
              mutate(Mean = mean_L3.F / Mean) %>%
              rename(Convert = Mean) %>%
              left_join(Fwds,conversion.rates.F, by = "Latest.Season") %>%
              mutate(Adj.Rate = ifelse(PP != 0,Rate * Convert,1/1000))
                  
  
      conversion.rates.D = Dmen %>%
                            group_by(Latest.Season) %>%
                            summarise(Mean = mean(Rate))
      
      mean_L3.D = conversion.rates.D %>%
                filter(Latest.Season == 2018) %>%
                .$Mean
      
      Dmen = conversion.rates.D %>%
              mutate(Mean = mean_L3.D / Mean) %>%
              rename(Convert = Mean) %>%
              left_join(Dmen,conversion.rates.D, by = "Latest.Season") %>%
              mutate(Adj.Rate = ifelse(PP != 0,Rate * Convert,1/1000))
      
  
  #Fitting priors to historical adjusted scoring rates

    #Forwards - fitting a Weibull prior
      
      #Plot Data
    
        h = hist(Fwds$Adj.Rate,
              col = "orchid4",
              border = "black",
              xlab = "Primary Points per Hour",
              main = "Fwds",
              breaks = 20)
      
        xfit = seq(min(Fwds$Adj.Rate), max(Fwds$Adj.Rate), length = 40)
    
      #Fit prior distribution
    
        F.Prior = fitdistr(Fwds$Adj.Rate,"weibull")
    
      #Overlay histogram with fitted Weibull distribution
    
        f.density = dweibull(xfit,shape = F.Prior$estimate[1], scale = F.Prior$estimate[2])
        f.density = f.density * diff(h$mids[1:2]) * length(Fwds$Adj.Rate)
        lines(xfit, f.density, col = "black", lwd = 2)
    
    #Defensemen - fitting a gamma prior
    
      #Plot Data
        
        h = hist(Dmen$Adj.Rate,
              col = "cadetblue4",
              border = "black",
              xlab = "Points per Hour",
              main = "Defensemen",
              breaks = 20)
        
        xfit = seq(min(Dmen$Adj.Rate), max(Dmen$Adj.Rate), length = 40)
    
      #Fit prior distribution
    
        D.Prior = fitdistr(Dmen$Adj.Rate,"gamma")
    
      #Overlay histogram with fitted gamma distribution
    
        d.density = dgamma(xfit,shape = D.Prior$estimate[1], rate = D.Prior$estimate[2])
        d.density = d.density * diff(h$mids[1:2]) * length(Dmen$Adj.Rate)
        lines(xfit, d.density, col = "black", lwd = 2)
    
  #FVF - A function for comparing posterior PP60 densities for forwards
  #Example... FVF("CONNOR.MCDAVID","SIDNEY.CROSBY")
    
  FVF = function(A,B){
    
    #Define Total Ice Time & Total Primary Points
      
        #Skater A
        
        TOI.A = (Data_L3[Data_L3$Skater == A, c("TOI")]/60)
        PP.A = Data_L3[Data_L3$Skater == A, c("PP")]
      
        #Skater B
        
        TOI.B = (Data_L3[Data_L3$Skater == B, c("TOI")]/60)
        PP.B = Data_L3[Data_L3$Skater == B, c("PP")]
    
    #Compute Posterior to a Constant of Proportionality
      
        #Skater A
        
        posterior.A = function(q){
        (q^PP.A)*exp(-q*TOI.A)*(q^(F.Prior$estimate[1] - 1))*exp(-(q/F.Prior$estimate[2])^F.Prior$estimate[1])}
    
        #Skater B
        
        posterior.B = function(q){
        (q^PP.B)*exp(-q*TOI.B)*(q^(F.Prior$estimate[1] - 1))*exp(-(q/F.Prior$estimate[2])^F.Prior$estimate[1])}
    
    #Find Maximum Posterior Density
        
        #Skater A
        
        MaxRate.A = optimize(posterior.A,c(0,5),maximum = T)$maximum
        MaxDensity.A = posterior.A(MaxRate.A)
    
        #Skater B
        
        MaxRate.B = optimize(posterior.B,c(0,5),maximum = T)$maximum
        MaxDensity.B = posterior.B(MaxRate.B)
    
    #Draw 100,000 Samples for Each Skater w/ Rejection Method
    
        #Skater A
        
          #Generate 2 million uniform random variables (min = 0, max = 5)
          
          Draws.A = as.data.frame(c(runif(2000000,min = 0,max = 5)))
          colnames(Draws.A)[1] = "Samples"
        
          #Evaluate posterior at each value generated by uniform distro  
          
          Draws.A$Eval = posterior.A(Draws.A$Samples)
          
          #Define rejection region  
          
          Draws.A$Rej.Region = Draws.A$Eval / MaxDensity.A
        
          #Generate 2 million uniform random variables (min = 0, max = 1)
          #and reject those which fall within the rejection region
          
          Draws.A$Rej.Test = c(runif(2000000,min =0,max = 1))
          Draws.A = Draws.A %>%
                      filter(Rej.Test < Rej.Region)
        
          #Keep first 100,000 posterior samples
          
          Draws.A = Draws.A[1:100000,]
    
        #Skater B
        
          #Generate 2 million uniform random variables (min = 0, max = 5)
          
          Draws.B = as.data.frame(c(runif(2000000,min = 0,max = 5)))
          colnames(Draws.B)[1] = "Samples"
          
          #Evaluate posterior at each value generated by uniform distro 
          
          Draws.B$Eval = posterior.B(Draws.B$Samples)
          
          #Define rejection region
          
          Draws.B$Rej.Region = Draws.B$Eval / MaxDensity.B
          
          #Generate 2 million uniform random variables (min = 0, max = 1)
          #and reject those which fall within the rejection region
          
          Draws.B$Rej.Test = c(runif(2000000,min =0,max = 1))
          Draws.B = Draws.B %>%
                      filter(Rej.Test < Rej.Region)
          
          #Keep first 100,000 posterior samples
          
          Draws.B = Draws.B[1:100000,]
          
    #Combine Skater Samples
        
        Draws = cbind.data.frame(Draws.A$Samples,Draws.B$Samples)
        
    #Estimate Probability A > B
        
        AvsB = round(nrow(Draws[Draws[,1] >= Draws[,2],]) / 1000,0)
    
    #Plot Posterior Distributions
        
        data = data.frame(dens = c(Draws.A$Samples, Draws.B$Samples)
                      ,Skaters = rep(c(A,B), each = 100000))
         
        plot.title = paste("5on5 Scoring Talent (Last 3 Seasons)"," ","\nProbability ",A," > ",B," = ",AvsB,"%",sep="")

        ggplot(data, aes(x = dens, fill = Skaters)) + 
          geom_density(alpha = 0.5) +
          ggtitle(plot.title) +
          theme(plot.title = element_text(hjust = 0.5), legend.pos = "bottom") +
          labs(x = "Scoring Rate (PrimaryP/60)", y = "Relative Likelihood")
  }
  
  #DVD - A function for comparing posterior PP60 densities for defensemen
  #Example... DVD("ROMAN.POLAK","RASMUS.RISTOLAINEN")
  
  DVD = function(A,B){
    
    #Define Total Ice Time & Total Primary Points
    
        #Skater A
    
        TOI.A = (Data_L3[Data_L3$Skater == A, c("TOI")]/60)
        PP.A = Data_L3[Data_L3$Skater == A, c("P")]
        
        #Skater B
      
        TOI.B = (Data_L3[Data_L3$Skater == B, c("TOI")]/60)
        PP.B = Data_L3[Data_L3$Skater == B, c("P")]
    
    #Draw 100,000 Samples for Each Skater  
        
        #Skater A
        
        Draws.A = as.data.frame(c(rgamma(100000,shape = D.Prior$estimate[1] + PP.A, rate = D.Prior$estimate[2] + TOI.A)))
        colnames(Draws.A)[1] = "Samples"
    
        #Skater B
        
        Draws.B = as.data.frame(c(rgamma(100000,shape = D.Prior$estimate[1] + PP.B, rate = D.Prior$estimate[2] + TOI.B)))
        colnames(Draws.B)[1] = "Samples"
        
    #Estimate Probability A > B
        
        Draws = cbind.data.frame(Draws.A$Samples,Draws.B$Samples)
        AvsB = round(nrow(Draws[Draws$`Draws.A$Samples` >= Draws.B$Samples,]) / 1000,0)
    
    #Plot Posterior Distributions
        
        data = data.frame(dens = c(Draws.A$Samples, Draws.B$Samples)
                       , Skaters = rep(c(A,B), each = 100000))
    
        plot.title = paste("5on5 Scoring Talent (Last 3 Seasons)"," ","\nProbability ",A," > ",B," = ",AvsB,"%",sep="")
    
        ggplot(data, aes(x = dens, fill = Skaters)) + 
          geom_density(alpha = 0.5) +
          ggtitle(plot.title) +
          theme(plot.title = element_text(hjust = 0.5), legend.pos = "bottom") +
          labs(x = "Scoring Rate (PrimaryP/60)", y = "Relative Likelihood")
  }
  
  
  #ggsave(filename = "DEFAULT.png",plot = last_plot(),dpi = 300, path = "C:/Users/mimic/Desktop/R plots")