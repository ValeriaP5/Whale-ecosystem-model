


title: "Stage structured models"

runtime: shiny
output: html_document
---

April 2021

```{r setup, include=FALSE}
library(shiny)
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
enableBookmarking(store = "url")

varnames <- rbind(c("Transition ages","age_trans"),
                          c("Natural survival", "p_s"),
                          c("Harvest rate","h_s"),
                          c("Fecundities", "births"))
```

```{r}
library(ggplot2)
use_plotly <- require(plotly)

survparams <- function(nstages, age_transition, p_stage, h_stage){
  
  x <- (nstages != length(age_transition)) | (nstages != length(p_stage)) |
    (nstages != length(h_stage)) 
  
  if(any(x)) stop("Number of stages and lengths of demographic parameters don't match")
  
  P <- rep(NA, nstages)
  G <- rep(NA, nstages)
  
  for (i in 1:nstages){
     
     p_natural_instant <- log(p_stage[i]) + log(1-h_stage[i])
     p_natural <- exp(p_natural_instant)
     if (i>1){
      years_in_stage <- age_transition[i] - age_transition[i-1]
     } else {
       years_in_stage <- age_transition[i]
     }
     
     P[i] <- ((1-p_natural^(years_in_stage-1))/(1-p_natural^(years_in_stage)))*p_natural
     G[i] <- ((p_natural^(years_in_stage)) * (1-p_natural))/(1-p_natural^(years_in_stage))
  }
 return(data.frame(P = P, G = G))
}

createA <- function(nstages, births, age_trans, p_s, h_s){
    x <- survparams(nstages, age_trans, p_s, h_s)
    A <- diag(nstages)*0
    A[1,] <- births
    diag(A) <- x$P
    irow <- 2:nstages
    icol <- 1:(nstages-1)
    rowcolind <- cbind(irow, icol)
    A[rowcolind] <- x$G[1:(nstages - 1)]
    colnames(A) <- paste("Age ", 1:nstages)
    A
}

 createX <- function(input, A){
    xtemp <- matrix(0, 
                    nrow = input$nstages, ncol = 3)
    xtemp[,1] <- asnum(input$stage_sizes)  
    xtemp
  }

asnum <- function(x){
  as.numeric(unlist(strsplit(x, ",")))
}

leslie <- function(tmax, Ninit, A){
  if (length(Ninit) != nrow(A)) stop("Dimension mismatch between initial abundances and Leslie matrix")
  
  N <- matrix(0, nrow = length(Ninit), ncol = tmax)
  N[,1] <- Ninit
  for (t in 1:(tmax-1)){
    N[,t+1] <- A %*% N[,t]
  }
  N
  
}
```

### Leslie matrix   

```{r}
inputPanel(
  numericInput("nstages", label = "Number of stages",
               value = 3, 
              min = 1, max = 20, step = 1),
  textInput("age_trans", label = "Transition ages", 
            value = "1, 15, 50"),
  textInput("p_s", label = "Natural survival probablities", 
            value = "0.5, 0.73, 0.97"),
  textInput("h_s", label = "Harvest mortality rates",
               value = c("0,0,0")), 
  textInput("births", label = "Fecundities", 
            value = "0,0,2"), 
  textInput("stage_sizes", label = "Stage abundances", 
            value = "1500, 4890, 750"),
  numericInput("tmax", label = "Number of years",
               value = 50, 
              min = 1, max = 200, step = 1),
  numericInput("yr_init", label = "Initial year",
               value = 2017, 
              min = 0, max = 2100, step = 1),
   checkboxInput("logged", label = "Log y-axis?", 
                 value = FALSE),
   actionButton("create_leslie", label = "Run simulation"),
    bookmarkButton()
)
```
```{r}

ptable <- eventReactive(input$create_leslie, {
  A <- createA(input$nstages, asnum(input$births), asnum(input$age_trans), 
               asnum(input$p_s), asnum(input$h_s))
  A
  })

xtable <- eventReactive(input$create_leslie, {
   createX(input, ptable())
  })


renderTable({
  ptable()
}, digits = 3)

renderText({
  paste0("The finite rate of growth for your Leslie matrix is ",
        signif(eigen(ptable())$values[1],4), ". \n", 
        "Estimated as the leading eigenvalue of the 
        Leslie matrix. ")
})
```

### Initial population structure 

```{r}

renderTable({
  xtemp <- xtable()
  xtemp[,2] <- ptable() %*% xtemp[,1]
  xtemp[,3] <- ptable() %*% xtemp[,2]
  
  colnames(xtemp) <- paste("Year", c(1,2,3))
  xtemp
})
```

### Temporal simulations

```{r}
N1 <- eventReactive(input$create_leslie, {
  leslie(input$tmax, xtable()[,1], ptable())
})

p1 <- eventReactive(input$create_leslie, {
  
  dplot <- data.frame(years = (1:input$tmax) + input$yr_init - 1,
                      Ntot = colSums(N1())
  )
  #Make long format data 
  Ndat <- data.frame(years = rep(dplot$years, each = input$nstages), 
                     stage = factor(rep(1:input$nstages, input$tmax)),
                     N = as.numeric(N1()))

   g1 <- ggplot(dplot, aes(x = years, y = Ntot)) + 
     geom_line() +
     geom_point() +
     ylab("Abundance") + 
     theme_classic()
   
   g2 <- ggplot(Ndat, aes(x = years, y = N, color = stage)) + 
     geom_line() + 
     geom_point() +
    ylab("Abundance") + 
     theme_classic() + 
     scale_color_brewer(palette = "Dark2")
   
   # cat(file=stderr(), "\n",which(Ndat$years %in% baryears), "stopped", "\n")
   
   if (input$logged){
     g1 <- g1 + scale_y_log10()
     g2 <- g2 + scale_y_log10()
   }
   
   return(list(g1 = g1, g2= g2))

   
})

```

#### Total abundance
```{r}
 if (use_plotly){
   
  renderPlotly({
    ggplotly(p1()[[1]])
  })
   
 } else {
   renderPlot({
     p1()[[1]]
   })
 }

```

#### Abundance by stage 


```{r}
 if (use_plotly){
   
  renderPlotly({
    ggplotly(p1()[[2]])
  })
   
 } else {
   renderPlot({
     p1()[[2]]
   })
 }

```

## Sensitivity analysis - temporal plot 

Select one parameters to vary and which lifestage the parameter should be varied for. For the first parameter enter minimum and maximum values, separated by a comma. The app will plot the total population size, with the different parameter values as different colours. 

```{r}
inputPanel(
  selectInput("param1", label = "First parameter to vary", 
              choices = c("Transition ages" = "age_trans",
                          "Natural survival" = "p_s",
                          "Harvest rate" = "h_s",
                          "Fecundities" = "births"), 
              selected = "p_s"),
  numericInput("param1_stage", label = "Parameter 1 life stage",
               value = 1),
  numericInput("param1_n", label = "Number of values to simulate",
               value = 5),
  textInput("param1_values", label = "Parameter 1 values",
            value ="0.5, 0.9"),
   actionButton("run_model_sens1", label = "Run simulations")
)
```

```{r}

psens <- eventReactive(input$run_model_sens1, {
  nstages <- input$nstages
  births <- asnum(input$births)
  age_trans <- asnum(input$age_trans)
  p_s <- asnum(input$p_s)
  h_s <- asnum(input$h_s)
  param1range <- asnum(input$param1_values)
  
  param1vals <- seq(param1range[1], param1range[2], 
                length.out = input$param1_n)
  
  Ndat <- NULL
  
  for (i in 1:input$param1_n){
    
    if (input$param1 == "p_s") p_s[input$param1_stage] <- param1vals[i]
      
    if (input$param1 == "h_s") h_s[input$param1_stage] <- param1vals[i]
    
    if (input$param1 == "age_trans") age_trans[input$param1_stage] <- round(param1vals[i])
      
    if (input$param1 == "births") births[input$param1_stage] <- param1vals[i]
      
      Atemp <- createA(nstages, births, age_trans, 
               p_s, h_s)
      Xtemp <- createX(input, ptable())
      Ntemp <- leslie(input$tmax, Xtemp[,1], Atemp)

  dplot <- data.frame(years = (1:input$tmax) + input$yr_init - 1,
                      Ntot = colSums(Ntemp),
                      Parameter_value = param1vals[i])

  Ndat <- c(Ndat, list(dplot))
  
  }
  
  Nout <- do.call("rbind", Ndat)
  var1 <- varnames[match(input$param1, varnames[,2]),1]
  g1 <- ggplot(Nout, aes(x = years, y = Ntot, color = Parameter_value)) + 
     geom_line() + 
     geom_point() +
    ylab("Abundance") + 
    xlab("Year") + 
     theme_classic() +
    labs(color = var1)
  return(g1)
})

```


```{r}
if (use_plotly){
  
  renderPlotly({
    ggplotly(psens())
  })
  
} else {
  renderPlot({
    psens()
  })
}

```  


## Sensitivity analysis - growth rates 

Select two parameters to vary and which lifestage the parameters should be varied for. For the first parameter enter minimum and maximum values, separated by a comma. For the second parameter, enter as many discrete values as you wish, seperated by commas. The app will plot the finite growth rate, with the first parameter on the x-axis and the second parameter as different colours. 

```{r}
inputPanel(
  selectInput("param2", label = "First parameter to vary", 
              choices = c("Transition ages" = "age_trans",
                          "Natural survival" = "p_s",
                          "Harvest rate" = "h_s",
                          "Fecundities" = "births"), 
              selected = "p_s"),
  numericInput("param2_stage", label = "Parameter 1 life stage",
               value = 1),
  numericInput("param2_n", label = "Number of values to simulate",
               value = 30),
  textInput("param2_values", label = "Parameter 1 values",
            value ="0.5, 0.9"),
  selectInput("param3", label = "Second parameter to vary",
              choices = c("Transition ages" = "age_trans",
                          "Natural survival" = "p_s",
                          "Harvest rate" = "h_s",
                          "Fecundities" = "births"), 
              selected = "births"),
  numericInput("param3_stage", label = "Parameter 2 life stage",
               value = 3),
  textInput("param3_values", label = "Parameter 2 values", 
            value = "2, 4, 6"),
   actionButton("run_model_sens2", label = "Run simulations")
)
```

```{r}
psens2 <- eventReactive(input$run_model_sens2, {
  nstages <- input$nstages
  births <- asnum(input$births)
  age_trans <- asnum(input$age_trans)
  p_s <- asnum(input$p_s)
  h_s <- asnum(input$h_s)
  param2range <- asnum(input$param2_values)
  
  param2vals <- seq(param2range[1], param2range[2], 
                length.out = input$param2_n)
  # cat(file=stderr(), "x", param1vals, "\n")
  param3vals <- asnum(input$param3_values)
  nparam3 <- length(param3vals)
  
  datparams <- NULL
  
  for (i in 1:input$param2_n){
    
    if (input$param2 == "p_s") p_s[input$param2_stage] <- param2vals[i]
      
    if (input$param2 == "h_s") h_s[input$param2_stage] <- param2vals[i]
    
    if (input$param2 == "age_trans") age_trans[input$param2_stage] <- round(param2vals[i])
      
    if (input$param2 == "births") births[input$param2_stage] <- param2vals[i]
      
    for (j in 1:nparam3){
      if (input$param3 == "p_s") p_s[input$param3_stage] <- param3vals[j]
      
      if (input$param3 == "h_s") h_s[input$param3_stage] <- param3vals[j]
      
      if (input$param3 == "age_trans") age_trans[input$param3_stage] <- round(param3vals[j])
      
      if (input$param3 == "births") births[input$param3_stage] <- param3vals[j]
      
      
      A <- createA(nstages, births, age_trans, 
                   p_s, h_s)
      
      dattemp <- data.frame(param2 = param2vals[i],
                            param3 = param3vals[j],
                            R = Re(eigen(A)$values[1]))
      datparams <- c(datparams, list(dattemp))
    }
  }
  
  dout <- do.call("rbind", datparams)
  

  var2 <- varnames[match(input$param2, varnames[,2]),1]
  var3 <- varnames[match(input$param3, varnames[,2]),1]
  g1 <- ggplot(dout, aes(x = param2, y = R, color = param3)) + 
    geom_point() +
    ylab("Finite growth rate per year") + 
    xlab(var2) + 
    theme_classic() +
    labs(color = var3) + 
    geom_hline(yintercept = 1)
  
  return(g1)
})

```


```{r}
if (use_plotly){
  
  renderPlotly({
    ggplotly(psens2())
  })
  
} else {
  renderPlot({
    psens2()
  })
}

```  
