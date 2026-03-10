#--------------------------------------#
# R shiny Reference Box
#--------------------------------------#

referenceContent <- reactive({
  switch(input$ref,
         "Ransick et al. (mouse, recommended)" = {
           content <- "<em>Ransick et al., 2019:</em><br>
              <p>Murine kidney cell atlas generated with Illumina Hi-Seq sequencing by 10X Genomics Chromium platform.
              Kidneys were derived from 2 adult male and 2 adult female C57BL6/J mice. 
              High resolution atlas as prior to cell dissociation kidney was subdivided into cortex, the outer medulla and the inner medulla. 
              30 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types.</p><br>
              <p>Data was downloaded from <a>https://github.com/qinzhu/kidneycellexplorer/tree/master/data</a> (2023/02/03). 
              Detailed description of the cell types displayed in the legend was downloaded from: <a>https://cello.shinyapps.io/kidneycellexplorer/</a> (2023/02/20)</p>
              <br><p>Reference: <em> Ransick, A., Lindström, N.O., Liu, J., Zhu, Q., Guo, J.-J., Alvarado, G.F., Kim, A.D., Black, H.G., Kim, J., McMahon, A.P., 2019. 
              Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. Dev. Cell 51, 399-413.e7. <a>https://doi.org/10.1016/j.devcel.2019.10.005</a>
              </em></p>"
         },
         "Park et al. (mouse)" = {
           content <- "<em>Park et al., 2018:</em><br>
              <p>First single cell RNA sequencing atlas of the mouse kidney generated with droplet-based single cell RNA sequencing. 
              Kidneys were derived from 7 healthy male mice. 
              24 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types.</p><br>
              <p>Data was downloaded from Gene Expression omnibus (GEO; accession no. GSE107585) and further processed with Seurat v2.3.4 for normalization. Cell types <q>novel 1</q> and <q>novel 2</q> were removed.</p><br>
              <p>Reference: <em> Park, J., Shrestha, R., Qiu, C., Kondo, A., Huang, S., Werth, M., Li, M., Barasch, J., Suszták, K., 2018. 
              Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science 360, 758–763. <a>https://doi.org/10.1126/science.aar2131</a></em></p>"
         },
         "Lake et al. (human)" = {
           content <- "<em>Lake et al., 2023:</em><br>
              <p>Droplet-based single-cell RNA sequencing based on 10x Genomics Chromium platform with Illumina Hi-Seq sequencing. 28 kidney biopsies were derived from 26 healthy human donors.
              77 cell cluster were identified, including epithelial, endothelial, stromal, immune and neural cell types. </p><br>
              <p>The h5Seurat file of scRNA-seq data was downloaded from Kidney Cell Atlas website kpmp.org repository section on 2025/04/08. Cell type names and abbreviations were adapted from supplementary table 4 of <em>Lake et al.</em>.</p><br>
              <p>References: <br><em>https://www.kpmp.org0</em><br>
              <p><em>Lake et al., 2023. An atlas of healthy and injured cell states and niches 
              in the human kidney. Nature 619, 585–594.</em></p>"
         },
         "Zhang et al. (human)" = {
           content <- "<em>Zhang et al., 2021:</em><br>
              <p>Cell atlases generated from a cohort of renal clear cell carcinomas (RCCs) and benign adjacent kidney tissue of 9 patients using single-cell RNA sequencing. 
              The benign scRNA-seq atlas demonstrated 16 cell clusters, including tubular epithelial, endothelial, stromal and immune cells. </p><br>
              <p>References: <br>
              <p><em>Zhang, Y. et al. Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response. Proceedings of the National Academy of Sciences 118, e2103240118 (2021).
                      </em></p>"
         }
  )
  content
})