library(shiny)


ui <- fluidPage(tabsetPanel(
  tabPanel(title = "GC content",
  plotOutput("GC")), 
  tabPanel(title = "Amino Acid Table",
           DT::dataTableOutput("AA")))
)

server <- function(input, output){
  output$GC <- renderPlot({
    library(seqinr)
    library(Biostrings)
    dengue <- read.fasta("./cov.fna")
    dengueseq <- dengue[[1]]
    
    starts <- seq(1, length(dengueseq)-2000, by = 2000)
    n <- length(starts)    # Find the length of the vector "starts"
    chunkGCs <- numeric(n) # Make a vector of the same length as vector
    
    for (i in 1:n) {
      chunk <- dengueseq[starts[i]:(starts[i]+1999)]
      chunkGC <- GC(chunk)
      chunkGCs[i] <- chunkGC
    }
    
    #"starts", but just containing zeroes
    title <- "Sliding window plot of GC content"

    plot(starts,chunkGCs,type="b",main = title, xlab="Nucleotide start position",ylab="GC content")
    })
  output$AA <- DT::renderDataTable({
    library(seqinr)
    library(Biostrings)
    dengue <- read.fasta("./cov.fna")
    dengueseq <- dengue[[1]]
    amin_acid <- translate(dengueseq)
    str <- NULL
    df <- read.table(text = "", colClasses = c("character", "integer"),
                     col.names = c("sequence", "length"))
    for (i in 1:length(amin_acid))
    {
      if ((amin_acid[i] == "*" || i == length(amin_acid)) && !is.null(str) )
      {
        df[nrow(df) + 1,] = list(str, nchar(str))
        str = NULL
      }
      else if (amin_acid[i] != "*")
        str <- paste(str, amin_acid[i], sep="")
    }
    df <- subset(df, !(df$length < 20))
    rownames(df) <- 1:nrow(df)
    df
  })
}

shinyApp(ui = ui, server = server)
