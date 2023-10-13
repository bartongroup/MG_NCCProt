# DIFFERENTIAL EXPRESSION EXPLORER
#
# Main UI
#

dashboardPage(
  
  dashboardHeader(
    title = "DE explorer"
  ),
  
  dashboardSidebar(
    div(
      style = {'padding: 10px'},
      h3(CONFIG$title),
      p(CONFIG$subtitle)
    ),
    mod_global_input_ui("global_input")
  ),
  
  dashboardBody(
    fluidRow(
      
      column(4,
             box(
               title = "Volcano/MA plot",
               status = "primary",
               width = NULL,
               height = NULL,
               mod_volma_plot_ui("volma_plot")
             ),
             box(
               title = "Feature plot",
               status = "primary",
               width = NULL,
               height = "auto",
               mod_feature_plot_ui("feature_plot")
             )
      ),
      
      column(4,
             box(
               title = "Feature information",
               status = "primary",
               width = NULL,
               height = 800,
               style = "overflow-y: scroll; overflow-x: hidden; height: 90%",
               mod_feature_info_ui("feature_info")
             )
      ),
      
      column(4,
             box(
               title = "Functional enrichment",
               status = "primary",
               width = NULL,
               height = 800,
               style = "overflow-y: scroll; overflow-x: hidden; height: 90%",
               mod_enrichment_ui("enrichment")
             )
      )
    )
  )
)
