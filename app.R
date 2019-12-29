library(MASS)
library(shiny)

# App for demonstration of the effect of various imputation methods on
# data with different missing data mechanisms: simulation using
# the bivariate normal distribution

ui <- fluidPage(
   titlePanel("Imputation of missing data - Demonstration"),
   
   sidebarLayout(
      sidebarPanel(
        # Allow LaTeX
        withMathJax(),
        # Distribution specifics
        numericInput("seed", 
                     "Random seed", 
                     72),
        numericInput("varx", 
                     "Variance of X",
                     1,
                     min = 10e-6,
                     max = 10e6),
        numericInput("vary", 
                     "Variance of Y",
                     1,
                     min = 10e-6,
                     max = 10e6),
        numericInput("cor", 
                     "Correlation",
                     0.5,
                     min = -1,
                     max = 1),
        
        # Choose missingness mechanism
        radioButtons("mechanism",
                    "Missingness mechanism",
                    c("MCAR", "MAR", "MNAR")),
        
        # MCAR: choose probability of observing a data point
        conditionalPanel(
         condition = "input.mechanism == 'MCAR'",
         sliderInput("percentage", 
                     "% of observed data",
                     0, 100, 75)),
        # MAR: define parameters that affect missingness
        conditionalPanel(
         condition = "input.mechanism == 'MAR'",
         numericInput("beta0", 
                      "\\( \\beta_0 \\)", 
                      1.5),
         numericInput("beta1", 
                      "\\( \\beta_1 \\)",
                      3)),
        # MNAR: define parameters that affect missingness
        conditionalPanel(
         condition = "input.mechanism == 'MNAR'",
         numericInput("beta0", 
                      "\\( \\beta_0 \\)", 
                      1.5),
         numericInput("beta1", 
                      "\\( \\beta_1 \\)",
                        3),
         numericInput("beta2", 
                      "\\( \\beta_2 \\)",
                      5)),
        
        # Choose analysis/imputation method
        selectInput("method",
                    "Method",
                    c("Complete case analysis", 
                      "Mean imputation", 
                      "Regression imputation",
                      "Stochastic regression imputation"))
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("imputed_plot"),
         #tableOutput('table')
      )
   )
)

server <- function(input, output) {
  # number of data points
  n <- 200

  # generate a bivariate normal distribution as specified and impose missingness in Y
  # according to pre-selected mechanism
  XY <- reactive({
    # use given random seed
    set.seed(input$seed)
    # covariance matrix
    Sigma <- matrix(c(input$varx, input$cor, input$cor, input$vary), nrow = 2, byrow = T)
    # draw a sample
    data <- mvrnorm(n = n, mu = c(0,0), Sigma = Sigma)
    
    # sample missingness indicator 
    if(input$mechanism == 'MCAR') {
      r <- rbinom(n, 1, input$percentage/100)
    }
    if(input$mechanism == 'MAR') {
      r <- rbinom(n, 1, 1/(1+exp(-(input$beta0+input$beta1*data[,1]))))
    }
    if(input$mechanism == 'MNAR') {
      r <- rbinom(n, 1, 1/(1+exp(-(input$beta0+input$beta1*data[,1]+input$beta2*data[,2]))))
    }
    
    # remove data where Y should be unobserved
    incomplete <- ifelse(r == 1, data[,2], NA)
    cbind(data[,1], incomplete)
  })
  
  # imputation according to the chosen method
  XY_imputed <- reactive({
    XY <- XY(); X <- XY[,1]; Y <- XY[,2]
    
    if(input$method == 'Complete case analysis') {
      XY
    } else {
      if(input$method == 'Mean imputation') {
        imputed <- ifelse(is.na(Y), mean(Y, na.rm=TRUE),Y)
      }
      if(input$method == 'Regression imputation') {
        reg <- lm(Y ~ X)
        imputed <- ifelse(is.na(Y), predict(reg, newdata = data.frame(X)), Y)
      }
      if(input$method == 'Stochastic regression imputation') {
        reg <- lm(Y ~ X)
        imputed <- ifelse(is.na(Y), predict(reg, newdata = data.frame(X)) + rnorm(n, 0, sigma(reg)), Y)
      }
      cbind(X, imputed)
    }
  }) 
   
  # visualisation
   output$imputed_plot <- renderPlot({
      plot(XY(), xlab="X", ylab="Y")
      points(XY_imputed()[is.na(XY()[,2]),], col = 2)
   })
   #output$table <- renderTable(table(...))
}

# Run the app 
shinyApp(ui = ui, server = server)

