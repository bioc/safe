`getLocalStatisticNames` <-
function(){
   return("safe.names" = c("f.ANOVA","t.LM","t.paired","t.Student","t.Welch","z.COXPH"),
          "gui.names"  = c("ANOVA F-statistic", "Regression",
                           "Paired Student's T-test", "Unpaired Student's T-test", 
                           "Welch T-test", "Cox Proportional Hazard Model"))
}

