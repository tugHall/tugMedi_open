###   I)  Define CLONE'S CLASSES ----------------------------------------------------------

#' Class 'Clone' for clones
#'
#' @keywords internal
clone <- setRefClass(

    Class = "Clone",

    fields = list(
        id = "numeric",          # identificator
        parent = "numeric",      # parent ID (for first - 0)
        N_cells = "numeric",     # Number of cells in clone
        c     = "numeric",       # Mean     of division counter
        c.var = "numeric",       # Variance of division counter
        d = "numeric",           # Probability of division
        i = "numeric",           # Probability of Hayflick limit
        a = "numeric",           # Probability of apoptosis
        k = "numeric",           # Probability of cell death by environment
        im = "numeric",          # Invasion/ metastasis probability; NaN for meta clones
        s = "numeric",           # Coefficient in the of apoptosis
        mutnum = "numeric",      # Number of poms             in Rother
        mutden = "numeric",      # Number of poms / exon size of Rother
        E = "numeric",           # Coefficient of friction term against to the split probability,
        Nmax = "numeric",        # The max number of cells that can exist in the primary tumor (Nmax = 1/E)

        Ha = "numeric",          # Apoptosis probability difference (apoptosis gene x weight)
        Him = "numeric",         # Invasion/ metastasis probability difference (invasion/ metastasis gene x weight)
        Hi = "numeric",          # Mitotic restriction probability difference (immortalized gene x weight)
        Hd = "numeric",          # Divide probability difference (cancer gene x weight)
        Hb = "numeric",          # Friction coefficient (angiogenic gene x weight)

        invasion     = "logical", # TRUE:           initiate invasion / metastasis
        invasion.suc = "numeric", # Num of cells to initiate invasion / metastasis
        primary  = "logical",    # Logical variable to show that primary transformation is already happened. 
                                 # True for primary and metastatic clones, FALSE for normal clones. 
        birthday = "numeric",    # Time step of birth of cell

        gene = "numeric",        # Flag for cancer gene function deficit (for drivers)
        pasgene = "numeric",     # Flag for cancer gene as passenger dysfunction
        PointMut_ID  =  "numeric", # ID of point mutation in object PointMut / pnt
        CNA_ID       =  "numeric", # ID of CNA mutation in object CNA

        drug.on = 'logical',
        drug.blocked.on = 'logical',
        drug.blocked.child = 'numeric', # id
        EF.mut_IDs = 'character',
        location.prmy = 'data.frame', # spatial info when primary    tumor
        location.meta = 'data.frame', # spatial info when metastatic tumor
        
        tmp = 'data.frame' # for anything temporal, should be overriden for re-use
    ),

    # Method
    methods = list(
        # Initialize
        initialize = function(gene_size, id=1, parent=0,
                              N_cells = 1,
                              c=0, c.var=0, d=dN, i=1,
                              a=0, k=kN, im=0,
                              s=sN, mutnum=0, mutden=0, E=1/K_N, Nmax=1/K_N,

                              invasion=FALSE, invasion.suc=0, primary=FALSE,
                              birthday=0,

                              gene=NULL, pasgene=NULL,
                              PointMut_ID = 0, CNA_ID = 0,

                              drug.on=FALSE, drug.blocked.on=FALSE, drug.blocked.child=-Inf,
                              EF.mut_IDs = '',
                              location.prmy = data.frame(x= 0, y= 0, z= 0, N_cells=NA),
                              location.meta = data.frame(x=NA, y=NA, z=NA, N_cells=NA) ) {
            id      <<- id
            parent  <<- parent
            N_cells <<- N_cells
            c     <<- c
            c.var <<- c.var
            d <<- d
            i <<- i
            if (is.null(a)) {
                a <<- .sigmoid( mutden * 1e6, s )
            } else {
                a <<- a
            }
            k  <<- k
            im <<- im
            s  <<- s
            mutnum <<- mutnum
            mutden <<- mutden
            E <<- E
            Nmax <<- 1.0 / E

            invasion     <<-  invasion
            invasion.suc <<-  invasion.suc
            primary      <<-  primary
            birthday     <<-  birthday

            Ha <<- 0
            Him <<- 0
            Hi <<- 0
            Hd <<- 0
            Hb <<- 0

            if (is.null(gene)) {
                gene <<- rep(0, gene_size)
            } else {
                gene <<- gene
            }
            if (is.null(pasgene)) {
                pasgene <<- rep(0, gene_size)
            } else {
                pasgene <<- pasgene
            }
            PointMut_ID <<-  PointMut_ID
            CNA_ID      <<-  CNA_ID

            drug.on            <<- drug.on
            drug.blocked.on    <<- drug.blocked.on
            drug.blocked.child <<- drug.blocked.child
            EF.mut_IDs <<- EF.mut_IDs
            location.prmy <<- location.prmy
            location.meta <<- location.meta
        },
        # Apoptosis
        calcApoptosis = function( onco1, Rother_genes, pnt_clones ) {
        
            .calcMutnum( clone1 = .self, onco1 = onco1, 
                         Rother_genes, pnt_clones )
        
            w = which( onco1$name %in% Rother_genes )
            if ( length( w ) > 0 ){
                L_CDS = sum( onco1$cds_1[ w ] + onco1$cds_2[ w ] ) / 2
                mutden <<- mutnum / L_CDS
                a      <<- .sigmoid( mutden * 1e6, s )
            } else {
                a <<- 0
            }
            
        }
    )
)


# 
.sigmoid <- function( xx, ss ) {

   yy <- 1 / ( 1 + exp( -1 * ss * ( xx - 0.5 ) ) )

   return( yy )
}


# Function to calculate mutnum as a number of poms in Rother region
.calcMutnum <- function( clone1, onco1, Rother_genes, pnt_clones ) {

    if ( length( Rother_genes ) == 0 ) {
        clone1$mutnum <- 0
    } else {
        mutnum <- 0 
        if ( length( pnt_clones ) != 0 && ( clone1$PointMut_ID != 0 )[ 1 ] ){
            for( id in clone1$PointMut_ID ){
                if ( pnt_clones[[ id ]]$Gene_name %in% Rother_genes ) mutnum = mutnum + 1
            }
        }
        clone1$mutnum <- mutnum
    }
}


#' Function to make one copy for clone1 in clone_init function
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return New object of class 'Clone' with the same info and new ID
#'
#' @keywords internal
clone_copy <- function( clone1 ) {
   pck.env$env$last_id = pck.env$env$last_id + 1

   clone2          <- clone1$copy()
   clone2$id       <- pck.env$env$last_id
   clone2$birthday <- pck.env$env$T
   clone2$parent   <- clone1$id

   return( clone2 )
}


#' Class 'Environ'
#'
#' @keywords internal
#'
environ <- setRefClass(
    # the class name
    Class = "Environ",

    # Fields
    fields = list(
        T = "numeric",           # time counter
        N = "numeric",           # number of normal cells
        P = "numeric",           # number of primary tumor cells
        M = "numeric",           # number of infiltrating / metastatic cells
        F = "numeric",           # a coefficient (Nmax = F/E) that determines the maximal number of cells
        c = "numeric",           # average number of divisions
        d = "numeric",           # mean value of splitting probability
        i = "numeric",           # average value of immortalization probability
        a = "numeric",           # average value of apoptosis probability
        k = "numeric",           # average probability of cell death by environment
        E = "numeric",           # average value of coefficients of friction term proportional to N, for splitting probability
        Nmax = "numeric",        # Maximal number of cells that can exist
        im = "numeric",          # average value of invasion / metastasis probability
        Ha = "numeric",          # average values of hallmarks:
        Him = "numeric",
        Hi = "numeric",
        Hd = "numeric",
        Hb = "numeric",
        type = "numeric",        # invasion / metastatic ratio
        gene = "numeric",        # cancer gene damage rate
        mutden = "numeric",      # density of mutations
        last_id = "numeric",     # last id of clone
        n_divisions = 'numeric'  # Number of divisions in a simulation
    ),

    # Methods
    methods = list(
        # Initialize
        initialize = function(Fb) {
            T <<- 0
            N <<- 0
            M <<- 0
            P <<- 0
            F <<- Fb
            n_divisions <<- 0
        }
    )
)


#' Class 'OncoGene'
#'
#' @keywords internal
oncogene <- setRefClass(
    Class = "OncoGene",

    fields = list(
        id = "numeric",         # id is same as in clone (key for clones)
        name  = "character",    # Cancer gene name list
        onsp  = "character",    # oncogene/suppressor indicator
        len   = "numeric",      # lengths of cancer genes
        cds_1 = "numeric",      # cancer gene CDS base length for parental chr 1
        cds_2 = "numeric",      # cancer gene CDS base length for parental chr 2
        rna_1 = "numeric",      # cancer gene RNA base number length for parental chr 1
        rna_2 = "numeric",      # cancer gene RNA base number length for parental chr 2
        p0_1  = "numeric",      # the probability of absent of mutations for parental chr 1
        p0_2  = "numeric",      # the probability of absent of mutations for parental chr 2
        prob_1 = "numeric",     # prob = c(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA) / sum(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA)
        prob_2 = "numeric",     # same for parental chr 2
        sum_prob_1 = "numeric",  # sum_prob = sum(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA)
        sum_prob_2 = "numeric"  # same for parental chr 2
    ),

    methods = list(
        read = function( Rint_file ) {

            name0 = NULL
            onsp0 = NULL
            cds0 = NULL
            rna0 = NULL

            data_int    =  read.table( Rint_file, sep = '\t', stringsAsFactors = FALSE, header = TRUE )
            if (  ! all( names( data_int ) ==  c( 'Gene', 'CDS_ID', 'Type' ) ) )
               stop( 'Column names in the Rint file should be Gene, CDS_ID, Type. ')

            for (i in 1:nrow(data_int)) {
                name <<- as.character(data_int[i, 'Gene'])
                if (!is.element(name, name0)) {
                    name0 = c(name0, name)
                    type = as.character(data_int[i, 'Type'])
                    if (type == "?") {
                        if (runif(1) > 0.5) {
                            type = "o"
                        } else {
                            type = "s"
                        }
                    }
                    onsp0 = c(onsp0, type)
                }
            }

            all_names  =  unique( pck.env$gene_map$Gene )
            add_names  =  all_names[ !( all_names  %in%  name0 ) ]

            name0  =  c( name0, add_names )
            for( gene in add_names ){
                if (runif(1) > 0.5) {
                    type = "o"
                } else {
                    type = "s"
                }
                onsp0 = c(onsp0, type)
            }

            # Get length of CDS and gene's length from gene_map:
            for ( i in 1:length(name0) ) {
                w     =  which( pck.env$gene_map$Gene == name0[ i ] )
                cds0  =  c( cds0, sum( pck.env$gene_map$End[w]  -  pck.env$gene_map$Start[w] + 1 ) )
                rna0  =  c( rna0, max( pck.env$gene_map$End[w] ) - min( pck.env$gene_map$Start[w]) + 1  )
            }

            id   <<- 1
            name <<- name0
            onsp <<- onsp0
            cds_1  <<- cds0
            cds_2  <<- cds0
            len  <<- length(name0)
            rna_1  <<- rna0
            rna_2  <<- rna0
            sum_prob_1 <<- sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
            sum_prob_2 <<- sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
            prob_1     <<- c( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob_1
            prob_2     <<- c( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob_2
            p0_1   <<- (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
            p0_2   <<- (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
        }
    )
)


#' Function to make one copy for onco1 in init_onco_clones function
#'
#' @param onco1 Object of class 'OncoGene'
#'
#' @return New object of class 'OncoGene' with the same info
#'
#' @keywords internal
onco_copy <- function( onco1 ){

   onco2 <- onco1$copy()

   return( onco2 )
}


#' Class 'Hallmarks'
#'
#' @keywords internal
#'
hallmark <- setRefClass(
    #
    Class = "Hallmarks",

    #
    fields = list(
        Ha = "numeric",       # Apoptosis hallmark indexes of genes in onco$name
        Hi = "numeric",       # Immortalization hallmark indexes of genes in onco$name.
        Hd = "numeric",       # Growth/antigrowth hallmark indexes of genes in onco$name.
        Hb = "numeric",       # Angiogenesis hallmark indexes of genes in onco$name.
        Him = "numeric",      # Invasion/metastatic transformation hallmark indexes of genes in onco$name.
        Ha_w = "numeric",     # Apoptosis hallmark weights of genes
        Hi_w = "numeric",     # Immortalization hallmark weights of genes
        Hd_w = "numeric",     # Growth/antigrowth hallmark weights of genes
        Hb_w = "numeric",     # Angiogenesis hallmark weights of genes
        Him_w = "numeric"     # Invasion/metastatic transformation hallmark weights of genes
    ),

    #
    methods = list(
        #
        read = function(file, normalization  =  TRUE, tumblers, compaction_factor ) {
            # normalization is an indicator to normalize all Hallmarks values
            data <- read.table( file, sep = "\t", header = TRUE, stringsAsFactors = FALSE )
            if ( ! all( names( data ) == c( 'Hallmark', 'Gene', 'Weight' ) ) )
               stop( 'Column names in the hallmark file should be Hallmark, Gene, Weight ' )

            Ha0 = NULL
            Hi0 = NULL
            HdN = NULL
            Hb0 = NULL
            Him0 = NULL
            Ha0_w = NULL
            Hi0_w = NULL
            HdN_w = NULL
            Hb0_w = NULL
            Him0_w = NULL

            # Acquire gene name and Hallmark coefficient by function from definition file
            for (i in 1:nrow(data)) {
                if (data[i, 'Hallmark' ] == "apoptosis") {
                    Ha0 = c(Ha0, as.character(data[i, 'Gene']))
                    Ha0_w = c(Ha0_w, as.numeric(as.character(data[i, 'Weight'])))

                } else if (data[ i, 'Hallmark' ] == "immortalization") {
                    Hi0 = c( Hi0, as.character(data[i, 'Gene' ] ) )
                    Hi0_w = c( Hi0_w, as.numeric( as.character( data[ i, 'Weight' ] ) ) )

                } else if ( data[ i, 'Hallmark' ] == "anti-growth" | data[i, 'Hallmark' ] == "growth" ) {
                    HdN = c( HdN, as.character( data[ i, 'Gene' ] ) )
                    HdN_w = c( HdN_w, as.numeric( as.character( data[ i, 'Weight' ] ) ) )

                } else if (data[i, 'Hallmark' ] == "angiogenesis") {
                    Hb0 = c( Hb0, as.character( data[ i, 'Gene' ] ) )
                    Hb0_w = c( Hb0_w, as.numeric( as.character( data[ i, 'Weight' ] ) ) )

                } else if (data[i, 'Hallmark' ] == "invasion") {
                    Him0 = c( Him0, as.character( data[ i, 'Gene' ] ) )
                    Him0_w = c( Him0_w, as.numeric( as.character( data[ i, 'Weight' ] ) ) )

                }
            }

            names  =  pck.env$onco$name    # unique( data[ , 'Gene' ] )
            Ha <<- match(Ha0, names)
            Hi <<- match(Hi0, names)
            Hd <<- match(HdN, names)
            Hb <<- match(Hb0, names)
            Him <<- match(Him0, names)

            if ( length( Ha  ) == 0 ) Ha  <<- 0 
            if ( length( Hi  ) == 0 ) Hi  <<- 0 
            if ( length( Hd  ) == 0 ) Hd  <<- 0 
            if ( length( Hb  ) == 0 ) Hb  <<- 0 
            if ( length( Him ) == 0 ) Him <<- 0 
            
            if ( length( Ha0_w ) == 0 ) Ha0_w = 0 
            if ( length( Hi0_w ) == 0 ) Hi0_w = 0 
            if ( length( HdN_w ) == 0 ) HdN_w = 0 
            if ( length( Hb0_w ) == 0 ) Hb0_w = 0 
            if ( length( Him0_w ) == 0 ) Him0_w = 0 
            
            
            
            # Total by genetic mode
            if ( normalization ){
                Ha_sum =  ifelse( test = tumblers$apoptosis,       yes = sum(Ha0_w),  no = 1 )
                Hi_sum =  ifelse( test = tumblers$immortalization, yes = sum(Hi0_w),  no = 1 )
                Hd_sum =  sum(HdN_w)
                Hb_sum =  ifelse( test = tumblers$angiogenesis,    yes = sum(Hb0_w),  no = 1 )
                Him_sum = ifelse( test = tumblers$metastasis,      yes = sum(Him0_w), no = 1)
            } else {
                Ha_sum = 1
                Hi_sum = 1
                Hd_sum = 1
                Hb_sum = 1
                Him_sum = 1
            }

            Ha_w <<- ifelse( test = rep( tumblers$apoptosis, length( Ha0_w )),       
                             yes = Ha0_w / Ha_sum, 
                             no  = rep( 0, length( Ha0_w ) ) 
                             )
            
            Hi_w <<- ifelse( test = rep( tumblers$immortalization, length( Hi0_w )),    
                             yes = Hi0_w / Hi_sum, 
                             no  = rep( 0, length( Hi0_w ) )
                            )
            
            Hd_w <<- HdN_w / Hd_sum
            
            Hb_w <<- ifelse( test = rep( tumblers$angiogenesis, length( Hb0_w )),
                             yes = Hb0_w / Hb_sum, 
                             no  = rep( 0, length( Hb0_w ) )
                            )
            Him_w <<- ifelse( test = rep( tumblers$metastasis, length( Him0_w )),
                              yes = Him0_w / Him_sum, 
                              no  = rep( 0, length( Him0_w ) )
                            )

            if ( compaction_factor ){
                Ha_w  <<-  pck.env$CF$Ha  * Ha_w
                Hi_w  <<-  pck.env$CF$Hi  * Hi_w
                Hd_w  <<-  pck.env$CF$Hd  * Hd_w
                Hb_w  <<-  Hb_w                     #    pck.env$CF$Hb  * Hb_w
                Him_w <<-  pck.env$CF$Him * Him_w
            }

        },

        # Change the cell variables
        updateClone = function( clone1, onco1, Rother_genes, pnt_clones, tumblers ) {
            
            # Apoptosis
            if ( tumblers$apoptosis ){
                clone1$calcApoptosis( onco1, Rother_genes, pnt_clones )
                clone1$Ha  =  sum( clone1$gene[Ha] * Ha_w )
                clone1$a   =  clone1$a - clone1$Ha
                if ( clone1$a < 0 ) {
                    clone1$a = 0
                }
            } else {
                clone1$calcApoptosis( onco1, Rother_genes, pnt_clones )
                clone1$a   = 0
                clone1$Ha  = 0
            }
            
            # Immortalization
            if ( tumblers$immortalization ){
                clone1$Hi = sum( clone1$gene[Hi] * Hi_w )
                clone1$i = 1 - clone1$Hi
                if ( clone1$i < 0 ) {
                    clone1$i = 0
                }
            } else {
                clone1$i  = 0
                clone1$Hi = 0
            }
            
            # Angiogenesis
            if ( tumblers$angiogenesis ) { 
                clone1$Hb = sum( clone1$gene[Hb] * Hb_w )
            } else {
                clone1$Hb = 0
            }

            clone1$E     =   ifelse( ! tumblers$angiogenesis,
                                       1 / pck.env$K_N,
                                     ( 1 / pck.env$K_N ) * 10 ** ( - pck.env$Fb * clone1$Hb )
                                    )   
            clone1$Nmax  =  1.0 / clone1$E

            # Oncogene / suppressor, coupled with invasion / metastasis
            clone1$Hd  = sum( clone1$gene[Hd]  * Hd_w )
            
            if ( tumblers$metastasis ){
                clone1$Him = sum( clone1$gene[Him] * Him_w )
            } else {
                clone1$Him = 0
            }
            
            if ( ! clone1$invasion ) {
                clone1$d = dN + clone1$Hd * ( 1 - pck.env$env$P * clone1$E )
            } else {
                clone1$d = dN + clone1$Hd
            }
            if ( clone1$d < 0 ) { clone1$d = 0 }
            if ( clone1$d > 1 ) { clone1$d = 1 }

            if ( ! is.nan( clone1$im ) ) {
                clone1$im = clone1$Him
            }

        },
        # Change the environment variables
        updateEnviron = function(env, clones) {
            sum_cell(env, clones)
        }
    )
)


#' Function to update Hallmarks and variables after division or under initialization
#'
#' @keywords internal
update_Hallmarks <- function( clone1, onco1 ) {
    tumblers = list( apoptosis       = pck.env$tumbler_for_apoptosis_trial, 
                     angiogenesis    = pck.env$tumbler_for_angiogenesis_trial, 
                     immortalization = pck.env$tumbler_for_immortalization_trial, 
                     metastasis      = pck.env$tumbler_for_metastasis_trial )
    
    pck.env$hall$updateClone( clone1 = clone1, onco1 = onco1, 
                              Rother_genes = pck.env$Rother_genes,
                              pnt_clones   = pck.env$pnt_clones, 
                              tumblers = tumblers 
                              )
}


###  II) CNA and Point Mutations: CLASSES -------------------------------------------------

#' Class 'Point_Mutations'
#'
#' @keywords internal
#'
Point_Mutations <- setRefClass(
    #
    Class = "Point_Mutations",
    #
    fields = list(
        PointMut_ID    = "numeric",     # ID of point mutation
        Allele         = "character",   # A or B allele
        Parental_1or2  = "numeric",     # 1 or 2
        Chr            = "character",   # Chromosome name
        Ref_pos        = "numeric",     # Reference position
        Phys_pos       = "vector",      # Physical positions [from:to]
        Delta          = "vector",      # Delta of positions [from:to]
        Copy_number    = "numeric",     # Copy number of allele
        Gene_name      = "character",   # Gene's name
        MalfunctionedByPointMut = "logical",   # True or False
        mut_order      = "numeric",     # order of mutation to reproduce the gene_map data.frame
        Ovlp_CNA_ID    = "numeric",     # CNA_IDs which overlapped the pnt and changed Copy_number
        rst.ratio      = "numeric"      # Resist ratio against drug intervention
    ),

    #
    methods = list(
        # Initialization with default values
        initialize = function( ID = 1, ... ){
            PointMut_ID   <<- ID   
            Allele        <<- ""   
            Parental_1or2 <<-  0   
            Chr           <<-  ""  
            Ref_pos       <<-  0   
            Phys_pos      <<-  0   
            Delta         <<-  0   
            Copy_number   <<-  0   
            Gene_name     <<-  ""  
            MalfunctionedByPointMut <<- NA 
            mut_order     <<-  0 
            Ovlp_CNA_ID   <<-  0 
            rst.ratio     <<-  0
        },

        # Function to safe data to data frame df1
        safe = function() {
            df1 = data.frame( PointMut_ID = PointMut_ID,
                              Parental_1or2 = Parental_1or2,
                              Chr = Chr,
                              Ref_pos = Ref_pos,
                              Phys_pos = paste0( '[', paste(Phys_pos, collapse = ', ' ),  ']'),
                              Delta = paste0( '[', paste(Delta, collapse = ', ' ),  ']'),
                              Copy_number = Copy_number,
                              Gene_name = Gene_name,
                              MalfunctionedByPointMut = MalfunctionedByPointMut,
                              mut_order  =  mut_order,
                              Ovlp_CNA_ID  =  Ovlp_CNA_ID,
                              rst.ratio = rst.ratio,
                              stringsAsFactors = FALSE 
                             )

            return( as.data.frame( df1 ) )
        }

    )
)



#' Class 'CNA_Mutations'
#'
#' @keywords internal
CNA_Mutations <- setRefClass(
    Class = "CNA_Mutations",

    fields = list(
        CNA_ID	  = "numeric",             # ID of CNA mutation
        Parental_1or2	= "numeric",       # 1 or 2
        dupOrdel      = "character",       # dup or del
        Chr	          = "character",       # Chromosome name
        Ref_start	    = "numeric",       # Reference start position
        Ref_end	      = "numeric",         # Reference end position
        Gene_names	    = "character",     # Genes' name
        MalfunctionedByCNA  = "logical",   # True of False
        mut_order       = "numeric",       # Order of point mutations and CNAs
        rst.ratio       = "numeric"        # Resist ratio against drug intervention
    ),

    methods = list(
        # Initialization with default values
        initialize = function( ID = 1, ... ){
            CNA_ID        <<- ID   
            Parental_1or2 <<- 0    
            dupOrdel      <<- ""   
            Chr           <<- ""   
            Ref_start     <<- 0    
            Ref_end       <<- 0    
            Gene_names    <<- ""   
            MalfunctionedByCNA <<- NA 
            mut_order     <<- 0
            rst.ratio     <<- 0
        },

        # Function to safe data to data frame df1
        safe = function() {
            df1 = data.frame( CNA_ID = CNA_ID,
                              Parental_1or2 = Parental_1or2,
                              dupOrdel = dupOrdel,
                              Chr = Chr,
                              Ref_start = Ref_start,
                              Ref_end = Ref_end,
                              Gene_names = Gene_names,
                              MalfunctionedByCNA = MalfunctionedByCNA,
                              mut_order  =  mut_order,
                              rst.ratio = rst.ratio,
                              stringsAsFactors = FALSE 
                             )

            return( as.data.frame( df1 ) )
        }

    )
)


