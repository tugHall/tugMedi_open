

load_lib_functions4posts <- function( ){

  pkgs = list( 
     stringr = c( 'str_split', 'str_replace_all' ),
     futile.logger = c( 'flog.debug',
                        'flog.info',
                        'flog.warn',
                        'flog.error',
                        'flog.appender',
                        'flog.threshold'
                       )
              )

   for( pck in names( pkgs ) ){
      require( package = pck, character.only = TRUE, include.only = pkgs[[ pck ]])
   }

}


load_lib_functions4posts()
#flog.threshold( 'DEBUG' )
flog.debug( 'flog level in posts: DEBUG' )


