

## Color scheme for the ethnicities in UK Biobank
ethnicities <- c(
    'Other/Unknown',
    'African',
    'Any other Asian background',
    'Any other Black background',
    'Any other mixed background',
    'Any other white background',
    'Asian or Asian British',
    'Bangladeshi',
    'Black or Black British',
    'British',
    'Caribbean',
    'Chinese',
    'Do not know',
    'Indian',
    'Irish',
    'Mixed',
    'Other ethnic group',
    'Pakistani',
    'Prefer not to answer',
    'White',
    'White and Asian',
    'White and Black African',
    'White and Black Caribbean',
    'Not stated',
    'Asian or Asian British',
    'Black or Black British')

ethnicity2col <- c(
    '#969696', ## Other/Unknown
    '#08306B', ## African
    '#9E9AC8', ## Any other Asian background
    '#2171B5', ## Any other Black background
    '#CD8500', ## Any other mixed background  #<==== Clare changed this from #969696 to #CD8500 to distinguish from Unknown! <=== Changed back!
    '#74C476', ## Any other white background
    '#9E9AC8', ## Asian or Asian British
    '#6A51A3', ## Bangladeshi
    '#2171B5', ## Black or Black British
    '#00441B', ## British
    '#6BAED6', ## Caribbean
    '#D94801', ## Chinese
    '#969696', ## Do not know
    '#3F007D', ## Indian
    '#C8CE46', ## Irish
    '#969696', ## Mixed
    '#969696', ## Other ethnic group
    '#807DBA', ## Pakistani
    '#969696', ## Prefer not to answer
    '#74C476', ## White
    '#FFC822', ## White and Asian
    '#339999', ## White and Black African
    '#339999', ## White and Black Caribbean
    '#969696', ## Not stated
    '#9E9AC8', ## Asian or Asian British (same colour and shape as Any other Asian background)
    '#2171B5' ## Black or Black British (same colour and shape as Any other Black background)
    )
    
ethnicity2char <-  c(
    2, ## Other/Unknown
    6, ## African
    0, ## Any other Asian background
    3, ## Any other Black background
    4, ## Any other mixed background  #<==== Clare changed this from 2 to 4 to distinguish from Unknown!
    2, ## Any other white background
    0, ## Asian or Asian British
    3, ## Bangladeshi
    3, ## Black or Black British
    4, ## British
    1, ## Caribbean
    5, ## Chinese
    2, ## Do not know
    4, ## Indian
    1, ## Irish
    1, ## Mixed
    2, ## Other ethnic group
    5, ## Pakistani
    2, ## Prefer not to answer
    2, ## White
    3, ## White and Asian
    0, ## White and Black African
    2, ## White and Black Caribbean
    1, ## Not stated
    0, ## Asian or Asian British (same colour and shape as Any other Asian background)
    3  ## Black or Black British (same colour and shape as Any other Black background)
    )

names(ethnicity2col) <- ethnicities
names(ethnicity2char) <- ethnicities


transform.ethnicity <- function(Ethnicity) {
  Ethnicity[is.na(Ethnicity)] = "Other/Unknown"
  Ethnicity[Ethnicity==""] = "Other/Unknown"
  Ethnicity[Ethnicity=="."] = "Other/Unknown"
  Ethnicity[Ethnicity=="Do not know"] = "Other/Unknown"
  Ethnicity[Ethnicity=="Do_not_know"] = "Other/Unknown"
  Ethnicity[Ethnicity=="Prefer not to answer"] = "Other/Unknown"
  Ethnicity[Ethnicity=="Prefer_not_to_answer"] = "Other/Unknown"
  Ethnicity[Ethnicity=="Other ethnic group"] = "Other/Unknown"
  Ethnicity[Ethnicity=="Other_ethnic_group"] = "Other/Unknown"
  Ethnicity[Ethnicity=="Mixed"] = "Any other mixed background"
  Ethnicity[Ethnicity=="White"] = "Any other white background"    
  Ethnicity[Ethnicity=="Asian or Asian British"] = "Any other Asian background"
  Ethnicity[Ethnicity=="Asian_or_Asian_British"] = "Any other Asian background"
  Ethnicity[Ethnicity=="Black or Black British"] = "Any other Black background"
  Ethnicity[Ethnicity=="Black_or_Black_British"] = "Any other Black background"
  Ethnicity[Ethnicity=="Any_other_Asian_background"] = "Any other Asian background"
  Ethnicity[Ethnicity=="Any_other_Black_background"] = "Any other Black background"
  Ethnicity[Ethnicity=="Any_other_mixed_background"] = "Any other mixed background"
  Ethnicity[Ethnicity=="Any_other_white_background"] = "Any other white background"
  Ethnicity[Ethnicity=="White_and_Asian"] = "White and Asian"
  Ethnicity[Ethnicity=="White_and_Black_African"] = "White and Black African"
  Ethnicity[Ethnicity=="White_and_Black_Caribbean"] = "White and Black Caribbean"
  return(Ethnicity)
}


ethnicity2pop <- function(Ethnicity) {
  Ethnicity = as.character(Ethnicity)
  Ethnicity = transform.ethnicity(Ethnicity)
  return(Ethnicity)
}


# Categories as in the showcase.
transform.ethnicity2 <- function(Ethnicity) {
  Ethnicity[is.na(Ethnicity)] = "Not stated"
  Ethnicity[Ethnicity==""] = "Not stated"
  Ethnicity[Ethnicity=="."] = "Not stated"
  Ethnicity[Ethnicity=="Do not know"] = "Not stated"
  Ethnicity[Ethnicity=="Do_not_know"] = "Not stated"
  Ethnicity[Ethnicity=="Prefer not to answer"] = "Not stated"
  Ethnicity[Ethnicity=="Prefer_not_to_answer"] = "Not stated"
  Ethnicity[Ethnicity=="Other ethnic group"] = "Other ethnic group"
  Ethnicity[Ethnicity=="Other_ethnic_group"] = "Other ethnic group"
  Ethnicity[Ethnicity=="Mixed"] = "Mixed"
  Ethnicity[Ethnicity=="White"] = "White"    
  Ethnicity[Ethnicity=="Asian or Asian British"] = "Asian or Asian British"
  Ethnicity[Ethnicity=="Asian_or_Asian_British"] = "Asian or Asian British"
  Ethnicity[Ethnicity=="Black or Black British"] = "Black or Black British"
  Ethnicity[Ethnicity=="Black_or_Black_British"] = "Black or Black British"
  Ethnicity[Ethnicity=="Any_other_Asian_background"] = "Any other Asian background"
  Ethnicity[Ethnicity=="Any_other_Black_background"] = "Any other Black background"
  Ethnicity[Ethnicity=="Any_other_mixed_background"] = "Any other mixed background"
  Ethnicity[Ethnicity=="Any_other_white_background"] = "Any other white background"
  Ethnicity[Ethnicity=="White_and_Asian"] = "White and Asian"
  Ethnicity[Ethnicity=="White_and_Black_African"] = "White and Black African"
  Ethnicity[Ethnicity=="White_and_Black_Caribbean"] = "White and Black Caribbean"
  return(Ethnicity)
}


ethnicity2pop2 <- function(Ethnicity) {
  Ethnicity = as.character(Ethnicity)
  Ethnicity = transform.ethnicity2(Ethnicity)
  return(Ethnicity)
}



# based on the interim release ethnicty table.
ethnicityCats1 = list("White"=c("British","Irish","Any other white background"),
    "Asian"= c("Indian","Pakistani","Bangladeshi","Chinese","Any other Asian background"),
    "Black"=c("African","Caribbean","Any other Black background"),
    "Mixed"=c("White and Asian","White and Black African","White and Black Caribbean","Any other mixed background"),
    "Other/Unknown"=c("Other/Unknown"))

                                        # based on the question/answers in UKBiobank showcase: http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21000

ethnicityCats2 = list("White"=c("British","Irish","Any other white background"),
    "Asian" = c("Indian","Pakistani","Bangladeshi","Any other Asian background"),
    "Black" =c("African","Caribbean","Any other Black background"),
    "Mixed" =c("White and Asian","White and Black African","White and Black Caribbean","Any other mixed background"),
    "Chinese" =c("Chinese"),
    "Other/Unknown"=c("Other/Unknown"))


# includes Asian or Asian British categories and splits unknown. Use ethnicity2pop2()
ethnicityCats3 = list("White"=c("White","British","Irish","Any other white background"),
    "Asian or Asian British" = c("Asian or Asian British","Indian","Pakistani","Bangladeshi","Any other Asian background"),
    "Black or Black British" =c("Black or Black British","African","Caribbean","Any other Black background"),
    "Mixed" =c("Mixed","White and Asian","White and Black African","White and Black Caribbean","Any other mixed background"),
    "Chinese" =c("Chinese"),
    "Other/Unknown"=c("Other ethnic group","Not stated"))


getEthCat <- function(z,cats=ethnicityCats1){
    y = unlist(cats)
    u = unlist(sapply(y,function(i) {
        h=names(cats)[sapply(cats,function(x) i%in%x)]
        if(length(h)==0) return(NA) else return(h)
    }))
    names(u)=y
    k = u[z]
    return(k)
}

# getEthCat(c("British","Irish",NA))

## Genotype cluster colors

# old ones with blue green orange
geno.cols = c("#555753","#FF7F00","#4DAF4A","#377EB8")
geno.cols.lighter = c("#888A85","#FDB462","#B3DE69","#80B1D3")
# based on viridis colours
geno.cols=c( "#555753","#FDE725FF", "#21908CFF" ,"#440154FF")
# colourblind friendly according to http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png
geno.cols=c( "#555753", "#0072B2","#009E73" , "#E69F00") # reversed so that the reference allele (usually the major allele) is blue, not orange!
geno.cols=c( "#555753", "#E69F00", "#009E73" ,"#0072B2")

## function to add alpha
print("function!!!!!!")
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}


