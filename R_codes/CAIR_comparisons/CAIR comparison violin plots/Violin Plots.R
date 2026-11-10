## 0. Installing necessary packages ##
#############################################################################################
if (!require("plotly")) install.packages("plotly")
if (!require("brunnermunzel")) install.packages("brunnermunzel")
if (!require("data.table")) install.packages("data.table")
if (!require("processx")) install.packages("processx")

## 0.1. Loading packages
########################
library (plotly)
library (brunnermunzel)
library (data.table)
library (processx)

###########################
## 1. Preparing the data ##
#############################################################################################
cair_file <- "../Supplementary_materials/CAIR_supplementary_files/Complete_proteome_CAIRs.csv" # rev 26: the copy in the repository (run from a folder of the repository, e.g. Results)
df <- fread(cair_file, select = c(1,2,3,5))
colnames(df) <- c("org", "cair", "tax", "code")
# rev 26: the relative effects (RE), logit-scale CIs and BH-adjusted p-values are read from the
# output of "Bifurcation tests.R", so that the figure shows the same values as Table 2
bm_results <- fread("Bifurcation tests (proteome).csv")
each_plot_width = 70
two_digits <- function (x) formatC(x, format = 'g', digits = ifelse(x >= 0.995, 3, 2)) # rev 26: e.g. 0.996 is not shown as 1
assign('num_of_plots', 0)
y_spans <- list() # rev 26: vertical extent of each pair of violins, for the axes of each row (section 3.1)
soft_span <- function (v) { # extent drawn by plotly (spanmode 'soft': two Silverman bandwidths beyond the data)
  bw <- 1.059 * min(sd(v), IQR(v) / 1.349) * length(v)^(-1/5)
  c(min(v) - 2 * bw, max(v) + 2 * bw)
}

#######################################################################################
## 2. Creating a function to draw violin plots and attach Brunner Munzel test results##
#############################################################################################
violinplot <- function(test_num, grp1, grp2){
    num_of_plots <<- num_of_plots + 1
    y_spans[[num_of_plots]] <<- range(soft_span(df$cair[grp1]), soft_span(df$cair[grp2]))
    BMtest <- bm_results[bm_results$test == test_num, ] # rev 26 (was: brunnermunzel.test of the two groups)
    if (BMtest$p_adj < 0.0001){
      sig <- "<b>***</b>"
    } else if (BMtest$p_adj < 0.001){
      sig <- "<b>**</b>"
    } else if (BMtest$p_adj < 0.05){
      sig <- "<b>*</b>"
    } else {
      sig <- ""
    }
    
    test <- paste(test_num, sig, ' RE: ', # rev 26: relative effect (was "ES") and logit-scale CI
                  two_digits(BMtest$re),
                  '\nCI: ',
                  two_digits(BMtest$logit_low),
                  ' - ',
                  two_digits(BMtest$logit_high),
                  sep = '')
    
    p <- plot_ly(type = 'violin') %>%
      add_trace(y = df$cair[grp1],
                x0 = test,
                side = 'negative',
                width = each_plot_width,
                opacity = 0.85,
                pointpos = -0.4,
                points = 'suspectedoutliers',
                fillcolor = "#1587a0",
                box = list(visible = T, width = 0.2),
                meanline = list(visible = T, color="#f0f4fa"),
                line = list (width = 0.3, color = "#0b7791"),
                marker = list (size = 1.5, color="#0b7791",
                               line = list (outlierwidth = 0.15),
                               symbol="diamond")) %>% 
      
      add_trace(y = df$cair[grp2],
                x0 = test,
                side = 'positive',
                width = each_plot_width,
                opacity = 0.85,
                pointpos = 0.4,
                points = 'suspectedoutliers',
                fillcolor = "#52c0a3",
                box = list(visible = T, width = 0.2),
                meanline = list(visible = T, color="#f0f4fa"),
                line = list (width = 0.3, color = "#36806a"),
                marker = list (size = 1.5, color="#36806a",
                               line = list (outlierwidth = 0.15),
                               symbol="diamond")) %>% 
      
      layout(font = list(family = "Arial, sans-serif", size = 6.5),
            yaxis = list(zeroline = F, color = "#444444",
            tickfont = list(family = "Arial, sans-serif", size = 6.5)),
            showlegend = F)
}

################################
## 3. Runing comparison tests ##
#############################################################################################
# WARNING: The tests contating a group with less than 10
# observations are #NOT to be executed and thus #HASHED
#.......................................................
# The below codes are based on the tree of life
#.......................................................
#library(plyr)
#count(df, vars = "code") # to count organisms in each phylum
fig_rows = 10
fig <- subplot(
  violinplot('T1', df$code<=56, df$code>=57),
  violinplot('T2', df$code<=55, df$code==56),
  violinplot('T3', df$code<=51, df$code<=55 & df$code>=52),
  violinplot('T4', df$code<=42, df$code<=51 & df$code>=43),
  violinplot('T5', df$code<=34, df$code<=42 & df$code>=35),
  violinplot('T6', df$code<=18, df$code<=34 & df$code>=19),
  violinplot('T7', df$code<=10, df$code<=18 & df$code>=11),
  violinplot('T8', df$code<=6, df$code<=10 & df$code>=7),
  violinplot('T9', df$code<=4, df$code<=6 & df$code>=5),
  violinplot('T10', df$code==1, df$code<=4 & df$code>=2),
# violinplot('T11', df$code<=3 & df$code>=2, df$code==4),
# violinplot('T12', df$code==2, df$code==3),
# violinplot('T13', df$code==5, df$code==6),
  violinplot('T14', df$code==7, df$code<=10 & df$code>=8),
  violinplot('T15', df$code==8, df$code<=10 & df$code>=9),
  violinplot('T16', df$code==9, df$code==10),
  violinplot('T17', df$code==11, df$code<=18 & df$code>=12),
  violinplot('T18', df$code==12, df$code<=18 & df$code>=13),
  violinplot('T19', df$code==13, df$code<=18 & df$code>=14),
  violinplot('T20', df$code==14, df$code<=18 & df$code>=15),
  violinplot('T21', df$code<=16 & df$code>=15, df$code<=18 & df$code>=17),
# violinplot('T22', df$code==15, df$code==16),
# violinplot('T23', df$code==17, df$code==18),
# violinplot('T24', df$code==19, df$code<=34 & df$code>=20),
# violinplot('T25', df$code==20, df$code<=34 & df$code>=21),
  violinplot('T26', df$code<=22 & df$code>=21, df$code<=34 & df$code>=23),
# violinplot('T27', df$code==21, df$code==22),
  violinplot('T28', df$code<=32 & df$code>=23, df$code<=34 & df$code>=33),
  violinplot('T29', df$code==23, df$code<=32 & df$code>=24),
# violinplot('T30', df$code==24, df$code<=32 & df$code>=25),
# violinplot('T31', df$code==25, df$code<=32 & df$code>=26),
  violinplot('T32', df$code<=29 & df$code>=26, df$code<=32 & df$code>=30),
  violinplot('T33', df$code<=27 & df$code>=26, df$code<=29 & df$code>=28),
# violinplot('T34', df$code==26, df$code==27),
# violinplot('T35', df$code==28, df$code==29),
  violinplot('T36', df$code==30, df$code<=32 & df$code>=31),
# violinplot('T37', df$code==31, df$code==32),
  violinplot('T38', df$code==33, df$code==34),
  violinplot('T39', df$code<=36 & df$code>=35, df$code<=42 & df$code>=37),
  violinplot('T40', df$code==35, df$code==36),
  violinplot('T41', df$code==37, df$code<=42 & df$code>=38),
  violinplot('T42', df$code<=39 & df$code>=38, df$code<=42 & df$code>=40),
# violinplot('T43', df$code==38, df$code==39),
  violinplot('T44', df$code==40, df$code<=42 & df$code>=41),
# violinplot('T45', df$code==41, df$code==42),
  violinplot('T46', df$code<=49 & df$code>=43, df$code<=51 & df$code>=50),
  violinplot('T47', df$code<=48 & df$code>=43, df$code==49),
# violinplot('T48', df$code<=47 & df$code>=43, df$code==48),
  violinplot('T49', df$code<=45 & df$code>=43, df$code<=47 & df$code>=46),
# violinplot('T50', df$code<=44 & df$code>=43, df$code==45),
  violinplot('T51', df$code==43, df$code==44),
# violinplot('T52', df$code==46, df$code==47),
# violinplot('T53', df$code==50, df$code==51),
  violinplot('T54', df$code<=53 & df$code>=52, df$code<=55 & df$code>=54), # corrected in rev 26 (was "df$code>=42")
  violinplot('T55', df$code==52, df$code==53),
  violinplot('T56', df$code==54, df$code==55),
  violinplot('T57', df$code<=66 & df$code>=57, df$code>=67),
# violinplot('T58', df$code==57, df$code<=66 & df$code>=58),
  violinplot('T59', df$code<=60 & df$code>=58, df$code<=66 & df$code>=61),
  violinplot('T60', df$code==58, df$code<=60 & df$code>=59),
# violinplot('T61', df$code==59, df$code==60),
  violinplot('T62', df$code==61, df$code<=66 & df$code>=62),
# violinplot('T63', df$code==62, df$code<=66 & df$code>=63),
# violinplot('T64', df$code<=64 & df$code>=63, df$code<=66 & df$code>=65), # hashed in rev 26: 8 proteomes in the upper clade
# violinplot('T65', df$code==63, df$code==64),
  violinplot('T66', df$code==65, df$code==66),
  violinplot('T67', df$code==67, df$code>=68),
  violinplot('T68', df$code<=75 & df$code>=68, df$code<=92 & df$code>=76),
  violinplot('T69', df$code<=70 & df$code>=68, df$code<=75 & df$code>=71),
  violinplot('T70', df$code==68, df$code<=70 & df$code>=69),
# violinplot('T71', df$code==69, df$code==70),
# violinplot('T72', df$code==71, df$code<=75 & df$code>=72),
# violinplot('T73', df$code<=73 & df$code>=72, df$code<=75 & df$code>=74), # corrected (was "df$code>=62") and hashed in rev 26: 4 proteomes in the upper clade
# violinplot('T74', df$code==72, df$code==73),
  violinplot('T75', df$code==74, df$code==75), # corrected in rev 26 (was "df$code>=75")
  violinplot('T76', df$code==76, df$code>=77),
  violinplot('T77', df$code==77, df$code>=78),
# violinplot('T78', df$code==78, df$code>=79),
  violinplot('T79', df$code<=87 & df$code>=79, df$code<=92 & df$code>=88),
  violinplot('T80', df$code==79, df$code<=87 & df$code>=80),
  violinplot('T81', df$code<=86 & df$code>=80, df$code==87),
# violinplot('T82', df$code==80, df$code<=86 & df$code>=81),
  violinplot('T83', df$code<=82 & df$code>=81, df$code<=86 & df$code>=83),
# violinplot('T84', df$code==81, df$code==82),
# violinplot('T85', df$code==83, df$code<=86 & df$code>=84),
# violinplot('T86', df$code==84, df$code<=86 & df$code>=85),
# violinplot('T87', df$code==85, df$code==86),
  violinplot('T88', df$code<=90 & df$code>=88, df$code<=92 & df$code>=91),
# violinplot('T89', df$code<=89 & df$code>=88, df$code==90),
# violinplot('T90', df$code==88, df$code==89),
# violinplot('T91', df$code==91, df$code==92),
  shareY = T, nrows = fig_rows)

## 3.1. Vertical axes (rev 26: each row starts and ends at a labeled gridline; two decimals)
######################################################################################
row_of <- ceiling(seq_along(y_spans) / ceiling(num_of_plots / fig_rows))
y_axes <- list()
for (r in unique(row_of)) {
  span <- range(unlist(y_spans[row_of == r]))
  step <- ifelse(diff(span) > 0.09, 0.05, 0.02)
  y_axes[[ifelse(r == 1, 'yaxis', paste0('yaxis', r))]] <- list(range = round(c(step * floor(span[1] / step), step * ceiling(span[2] / step)), 2),
                                                                dtick = step, tick0 = 0, tickformat = '.2f')
}
fig <- do.call(layout, c(list(fig), y_axes))

##########################
## 4. Displaying output ##
#############################################################################################
# for more information on orca see: https://github.com/plotly/orca#Installation
orca(p = fig, file = "Fig2B.svg", format = 'svg', # Figure 2B of the article
     width = ceiling(num_of_plots/fig_rows)*each_plot_width, # rev 26: ceiling() instead of round(), so that 54 plots keep six columns
     height = fig_rows*each_plot_width)
