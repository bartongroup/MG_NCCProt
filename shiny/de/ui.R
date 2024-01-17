# DIFFERENTIAL EXPRESSION EXPLORER
#
# Main UI
#


page_sidebar(
  theme = bs_add_rules(
    bs_theme(bootswatch = "journal"),
    as_sass("
            div.dataTables_wrapper  div.dataTables_filter {
              width: 100%;
              float: none;
              text-align: left;
            }"
    )
  ),
  title = "DE explorer",

  sidebar = sidebar(
    title = CONFIG$title,
    # div(
    #   style = {'padding: 10px'},
    #   h3(CONFIG$title),
    #   p(CONFIG$subtitle)
    # ),
    mod_global_input_ui("global_input")
  ),
  
  layout_column_wrap(
    width = 1/3,
    layout_column_wrap(
      width = 1,
      card(mod_volma_plot_ui("volma_plot")),
      card(mod_feature_plot_ui("feature_plot"))
    ),
    card(mod_feature_info_ui("feature_info")),
    card(mod_enrichment_ui("enrichment")),
  )
)
  
