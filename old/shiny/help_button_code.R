#----------------------------------------------------#
# help buttons
#----------------------------------------------------#

observeEvent(input$help_test, {
  showModal(modalDialog(
    title = "How to upload your sample",
    HTML("To ensure accurate matching, make sure your datatable meets the following criteria:<br><br>
        <b>filetype:</b><br>.csv or .xlsx<br><br>
        <b>format:</b><br>
        <ul><li>Sample names in columns and Gene identifiers (Gene ID) in rows</li>
            <li>Column of Gene ID's is named: <b>Gene.names</b></li>
            <li>Table should not include any other columns except 'Gene.names' and sample data</li>
            <li>Select unit of bulk RNA-seq data: <b>counts</b> or <b>counts per million(CPM)</b></li>
            <li>If only transcripts per million (TPM) data is available select CPM button but be aware that matching results might be less accurate</li></ul>
        Please refer to the example below:"),
    img(src = "sample_datatable_example.png", height = 300, width = 425),
  ))
})

observeEvent(input$help_ref, {
  showModal(modalDialog(
    title = "References",
    HTML("We provide 2 murine single cell RNA-seq reference datasets and 1 human <b>single nucleus</b> RNA-seq reference dataset:<br><br>
        <ul><li>Ransick <em>et al.</em> from mouse (recommended)</li>
            <li>Park <em>et al.</em> from mouse </li>
            <li>Kidney Precision Medicine Project (KPMP) snRNA-seq from healthy human donors </li></ul><br>
        Ransick <em>et al.</em> consistently provided the most accurate results, followed by Park <em>et al.</em> and KPMP. <p>For more info on how we evaluated the methods, refer to the tab <q>About</q>.</p>
        For more info on how we obtained and prepared scRNA-seq references, refer to the tab <q>References</q>.")
  ))
})

observeEvent(input$help_gene, {
  showModal(modalDialog(
    title = "How to select geneset",
    "We provide 2 geneset options to run dataset comparisons:",
    br(),
    br(),
    strong("all genes:"),
    "compares gene expression profile between all genes contained in both reference and sample. ",
    br(),
    strong("marker genes tubular celltypes:"),
    "marker gene list containing 153 marker genes for tubular kidney epithelial cells derived from Lake", 
    em("et al."),
    "(1)",
    br(),
    br(),
    p("We recommend to start with global gene expression. For reference Ransick et al. and method Euclidean distance matching accuracy might improve when selecting tubular marker genes. For references Park et al. and KPMP, 
    global gene expression performed better."),
    br(),
    em("(1) Supplementary Table 5 from: Lake, B.B., et al., 2023. An atlas of healthy and 
    injured cell states and niches in the human kidney. Nature 619, 585–594")
  ))
})

observeEvent(input$help_upload, {
  showModal(modalDialog(
    title = "Upload your own marker gene list:",
    fluidRow(
      column(9,
             HTML("To upload your own set of marker genes refer to the example on the right:<br>
                  <ul><li>filetype: <b>.csv</b> or <b>.xlsx</b></li>
                      <li>file contains 1 column <q><b>Gene.names</b></q> with gene identifiers of your choice</li>
                  </ul><br>
                  Your marker gene list will subsequently be included in <q>3. Select your geneset</q>. You can upload as many files as you like
                  and then select the gene set of your choice.")),
      column(3,
             img(src = "marker_genes_format.png", height = 200, width = 110))
    )
  ))
})

observeEvent(input$help_demo, {
  showModal(modalDialog(
    title = NULL,
    HTML("To explore CellMatchR, we provide 3 preexisting bulk RNA-sequencing datasets for matching:<br><br>
         <ul><li><b>Nephron primary cells after kidney microdissection </b> by<em> Chen et al.</em></li>
         <li><b>2 replicates of human proximal tubule cell line HK-2</b> by<em> Khundmiri et al.</em></li>
         <li><b>4 replicates of mouse inner medullary collecting duct cell line m-IMCD3</b> by our collaborators Prof. Dr. Michael Köttgen and Dr. Lukas Westermann
         </ul><br>References:<br><em>Chen, L., Chou, C.-L., Knepper, M.A., 2021. A Comprehensive Map of mRNAs and Their Isoforms across All 14 Renal Tubule Segments of Mouse. J. Am. Soc. Nephrol. 32, 897. <a>https://doi.org/10.1681/ASN.2020101406</a></em>
         <br><br><em>Khundmiri, S.J., Chen, L., Lederer, E.D., Yang, C.-R., Knepper, M.A., 2021. Transcriptomes of Major Proximal Tubule Cell Culture Models. J. Am. Soc. Nephrol. 32, 86. <a>https://doi.org/10.1681/ASN.2020010009</a></em>"

  )))
})

observeEvent(input$help_heatmap, {
  showModal(modalDialog(
    title = "How the heatmap is created:",
    HTML("<ul><li>Sample and reference data are <b>log2 transformed</b></li>
          <li>Genes are clustered using <b>Ward D2</b> clustering method</li>
          <li>Cell types are <b>ordered along x-axis by descending rho</b> of Spearman's correlation</li>
          <li>If <q><b>scale by row</b></q> is selected, counts per million (CPM) of each gene (x) is scaled across all cell types using the following formula:</li>
          <div style ='text-align:center;'><b>x-min(x)/max(x)-min(x)</b></div>")
  ))
})

observeEvent(input$help_select, {
  showModal(modalDialog(
    title = "When to (un)select cell types:",
    HTML("<ul><li><b>Comparisons between reference and sample(s) are run between all selected cell types</b></li><li>For this reason, you have the option 
          to <b>unselect cell types</b> that you are not interested in both in reference and sample</li>
          <li>This is especially relevant if your uploaded datatable <b>contains RNA-seq data from different groups of samples</b> which you do not
          want to combine during matching</li></ul>")
  ))
})

observeEvent(input$rho_interpret, {
  showModal(modalDialog(
    title = NULL,
    HTML("<ul><li>Spearman's correlation compares the <b>ranks of gene expression</b> of the reference to the ranks of gene expression of the sample(s)</li>
          <li><b>The more similar the ranks, the higher the correlation coefficient rho</b></li>
          <li>Results show <b>Spearman's rho</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b></li>
          <li>Results are ordered from <b>highest to lowest rho</b>. <b>High rho means higher correlation</b> between reference and sample gene expression ranks</li>
          <li>High-dimensional genesets, i.e. all genes, lead to high rho's and small distinctions in between cell types</li>
          <li>Lower dimensional genesets, i.e. marker genes, reduce rho's and increase rho differences between cell types</li>
          <li>Testing has shown that in most cases, despite these differences, <b>the order of results</b> matters most to determine similarity between reference and sample(s)</li>
          <li><q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Spearman's rho of sample replicates</li>
          <li>The blue bars show the correlation in between the sample(s)</li></ul>")
  ))
})

observeEvent(input$ED_interpret, {
  showModal(modalDialog(
    title = NULL,
    HTML("<ul><li>Euclidean distance is a <b>geometric distance metric</b> that measures the distance of the reference genes to the sample genes in a 3D space</li>
          <li><b>Higher similarity of gene expression</b> of the reference to the sample <b>decreases the distance</b> of the datapoints of the genes in the geometric space</li>
          <li>Low Euclidean distance means higher similarity</b> between reference and sample gene expression</li> <li>Results show <b>Euclidean distance</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b></li>
          <li>Results are ordered from <b>lowest to highest Euclidean distance</b></li>
          <li>High-dimensional genesets, i.e. all genes, lead to lower distances and small distinctions in between cell types</li>
          <li>Lower dimensional genesets, i.e. marker genes, increase Euclidean distance and increase differences between cell types</li>
          <li>Testing has shown that in most cases, despite these differences, <b>the order of results</b> matters most to determine similarity between reference and sample(s)</li>
          <li>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Euclidean distance of all sample replicates</li>
         <li>The blue bars show the Euclidean distance between the sample(s)</li></ul></ul>")
  ))
})