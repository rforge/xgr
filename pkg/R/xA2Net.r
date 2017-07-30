#' Function to visualise an igraph object using ggnetwork
#'
#' \code{xA2Net} is supposed to visualise an igraph object using ggnetwork.
#'
#' @param g an object of class "igraph"
#' @param node.label either a vector labelling nodes or a character specifying which node attribute used for the labelling. If NULL (by default), no node labelling
#' @param label.wrap.width a positive integer specifying wrap width of node labelling
#' @param node.label.size a vector specifying node size or a character specifying which node attribute used for node label size
#' @param node.label.color the node label color
#' @param node.label.alpha the 0-1 value specifying transparency of node labelling
#' @param node.label.padding the padding around the labeled node
#' @param node.label.arrow the arrow pointing to the labeled node
#' @param node.label.force the repelling force between overlapping labels
#' @param node.shape an integer specifying node shape
#' @param node.xcoord a vector specifying x coordinates. If NULL, it will be created using igraph::layout_with_kk
#' @param node.ycoord a vector specifying y coordinates. If NULL, it will be created using igraph::layout_with_kk
#' @param node.color a character specifying which node attribute used for node coloring
#' @param node.color.title a character specifying the title for node coloring
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param ncolors the number of colors specified over the colormap
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param node.size either a vector specifying node size or a character specifying which node attribute used for the node size
#' @param node.size.title a character specifying the title for node sizing
#' @param node.size.range the range of actual node size
#' @param slim the minimum and maximum values for which sizes should be plotted
#' @param title a character specifying the title for the plot
#' @param edge.color a character specifying the edge colors
#' @param edge.color.alpha the 0-1 value specifying transparency of edge colors
#' @param edge.curve a numeric value specifying the edge curve. 0 for the straight line
#' @param edge.arrow.gap a gap between the arrow and the node
#' @return
#' a ggplot object
#' @note none
#' @export
#' @seealso \code{\link{xA2Net}}
#' @include xA2Net.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev/"
#' 
#' # load REACTOME
#' # restricted to Immune System ('R-HSA-168256') or Signal Transduction ('R-HSA-162582')
#' g <- xRDataLoader(RData.customised='ig.REACTOME', RData.location=RData.location)
#' neighs.out <- igraph::neighborhood(g, order=vcount(g), nodes="R-HSA-168256", mode="out")
#' nodeInduced <- V(g)[unique(unlist(neighs.out))]$name
#' ig <- igraph::induced.subgraph(g, vids=nodeInduced)
#'
#' # visualise the graph with vertices being color-coded
#' V(ig)$degree <- igraph::degree(ig)
#' gp <- xA2Net(g=ig, node.label='term_id', label.wrap.width=30, node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=1, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='degree', node.color.title='Degree', colormap='grey-orange-darkred', ncolors=64, zlim=c(0,10), node.size.range=3, edge.color="black",edge.color.alpha=0.3,edge.curve=0.05,edge.arrow.gap=0.02, title='')
#' gp <- xA2Net(g=ls_ig, node.label='term_id', label.wrap.width=30, node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=1, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='degree', node.color.title='Degree', colormap='grey-orange-darkred', ncolors=64, zlim=c(0,10), node.size.range=3, edge.color="black",edge.color.alpha=0.3,edge.curve=0.05,edge.arrow.gap=0.02, title='')
#' 
#' ###########################
#' # visualise gene network
#' glayout <- igraph::layout_with_kk(g)
#' V(g)$xcoord <- glayout[,1]
#' V(g)$ycoord <- glayout[,2]
#' V(g)$degree <- igraph::degree(g)
#' gp <- xA2Net(g=g, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0, node.label.arrow=0, node.label.force=0.01, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='priority', node.color.title='5-star\nrating', colormap='yellow-red', ncolors=64, zlim=c(0,5), node.size='degree', node.size.title='Degree', slim=c(0,5), edge.color="orange",edge.color.alpha=0.5,edge.curve=0,edge.arrow.gap=0.025, title='')
#' gp_rating <- xA2Net(g=g, node.label='name', node.label.size=2, node.label.color='black', node.label.alpha=0.8, node.label.padding=0.1, node.label.arrow=0, node.label.force=0.01, node.shape=19, node.xcoord='xcoord', node.ycoord='ycoord', node.color='priority', node.color.title='5-star\nrating', colormap='white-yellow-red', ncolors=64, zlim=c(0,5), node.size.range=5, edge.color="orange",edge.color.alpha=0.3,edge.curve=0,edge.arrow.gap=0.02, title='')
#' }

xA2Net <- function(g, node.label=NULL, label.wrap.width=NULL, node.label.size=NULL, node.label.color='darkblue', node.label.alpha=0.8, node.label.padding=1, node.label.arrow=0.01, node.label.force=1, node.shape=19, node.xcoord=NULL, node.ycoord=NULL, node.color=NULL, node.color.title=NULL, colormap='grey-orange-darkred', ncolors=64, zlim=NULL, node.size=NULL, node.size.title=NULL, node.size.range=c(1,4), slim=NULL, title='', edge.color="black", edge.color.alpha=0.5, edge.curve=0.1, edge.arrow.gap=0.02)
{
    
   	if(any(class(g) %in% c("igraph"))){
		ls_ig <- list(g)
	}else if(class(g)=="list"){
		## Remove null elements in a list
		ls_ig <- base::Filter(base::Negate(is.null), g)
		if(length(ls_ig)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'igraph' objects or a 'igraph' object.\n")
	}
	
	## check list_names
	ls_names <- names(ls_ig)
	if(is.null(ls_names)){
		ls_names <- paste('IG', 1:length(ls_ig), sep='_')
		names(ls_ig) <- ls_names
	}
    
	ls_df <- lapply(1:length(ls_ig), function(i){
    	ig <- ls_ig[[i]]
    	
		nnode <- igraph::vcount(ig)
		## node.xcoord (by default, NULL)
		if(length(node.xcoord)!=nnode | length(node.ycoord)!=nnode){
			if(!is.null(node.xcoord)){
				node.xcoord <- igraph::vertex_attr(ig, node.xcoord)
			}
			if(!is.null(node.ycoord)){
				node.ycoord <- igraph::vertex_attr(ig, node.ycoord)
			}
			
			if(is.null(node.xcoord) | is.null(node.ycoord)){
				## layout
				#glayout <- igraph::layout_with_kk(ig)
				glayout <- igraph::layout_as_tree(ig,root=dnet::dDAGroot(ig),circular=TRUE,flip.y=TRUE)
				glayout <- glayout[,c(2:1)]
				node.xcoord <- glayout[,1]
				node.ycoord <- glayout[,2]
			}
		}
		## scale into [-1,1]
		if(max(node.xcoord) != min(node.xcoord)){
			node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
		}
		if(max(node.ycoord) != min(node.ycoord)){
			node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
		}
	
		## node.label (by default, NULL)
		if(length(node.label)!=nnode){
			if(!is.null(node.label)){
				node.label <- igraph::vertex_attr(ig, node.label)
			}
			if(is.null(node.label)){
				node.label <- rep('', nnode)
			}
		}
		## text wrap
		if(!is.null(label.wrap.width)){
			width <- as.integer(label.wrap.width)
			res_list <- lapply(node.label, function(x){
				if(!is.na(x)){
					x <- gsub('_', ' ', x)
					y <- strwrap(x, width=width)
					if(length(y)==2){
						paste(y, collapse='\n')
					}else if(length(y)>2){
						#paste0(y[1], '...')
						paste0(paste(y[1:2],collapse='\n'),'...' )
					}else{
						y
					}
				}else{
					x
				}
			})
			node.label <- unlist(res_list)
		}
		V(ig)$n.label <- node.label
	
		## node.label.size (by default, 0)
		if(length(node.label.size)!=nnode){
			if(!is.null(node.label.size)){
				tmp.node.label.size <- igraph::vertex_attr(ig, node.label.size)
			}else{
				tmp.node.label.size <- rep(0, nnode)
			}
			if(is.null(tmp.node.label.size)){
				node.label.size <- rep(node.label.size, nnode)
			}else{
				node.label.size <- tmp.node.label.size
			}
		}
		V(ig)$n.label.size <- node.label.size
	
		## node.color (by default, 0)
		if(length(node.color)!=nnode){
			if(!is.null(node.color)){
				node.color <- igraph::vertex_attr(ig, node.color)
			}
			if(is.null(node.color)){
				node.color <- rep(0, nnode)
			}
		}
		## zlim
		if(is.null(zlim)){
			zlim <- c(min(node.color), max(node.color))
		}
		node.color[node.color<=zlim[1]] <- zlim[1]
		node.color[node.color>=zlim[2]] <- zlim[2]
		V(ig)$n.color <- node.color
	
		## node.size (by default, 0)
		if(length(node.size)!=nnode){
			if(!is.null(node.size)){
				tmp.node.size <- igraph::vertex_attr(ig, node.size)
			}else{
				tmp.node.size <- rep(0, nnode)
			}
			if(is.null(tmp.node.size)){
				node.size <- rep(node.size, nnode)
			}else{
				node.size <- tmp.node.size
			}
		}
		## slim
		if(is.null(slim)){
			slim <- c(min(node.size), max(node.size))
		}
		node.size[node.size<=slim[1]] <- slim[1]
		node.size[node.size>=slim[2]] <- slim[2]
		V(ig)$n.size <- node.size

		gnet <- ggnetwork::ggnetwork(intergraph::asNetwork(ig), layout=cbind(node.xcoord,node.ycoord), arrow.gap=edge.arrow.gap, cell.jitter=0.75)
		data.frame(gnet, trait=rep(names(ls_ig)[i],nrow(gnet)), stringsAsFactors=F)
	})
    df <- do.call(rbind, ls_df)
    
    ## To replace only factors:
    i <- sapply(df, is.factor)
	df[i] <- lapply(df[i], as.character)

    ## ordered according to the input
    df$trait <- factor(df$trait, levels=names(ls_ig))

    ## make sure numeric
    df$n.color <- as.numeric(df$n.color)
    df$n.size <- as.numeric(df$n.size)
    
    #############################################################
    n.color <- n.size <- n.label <- n.label.size <- NULL
    x <- y <- xend <- yend <- NULL
    
	## ggplot
	gp <- ggplot(df, aes(x=x,y=y,xend=xend,yend=yend)) 
	
	if(igraph::is_directed(ls_ig[[1]])){
		gp <- gp + ggnetwork::geom_edges(color=edge.color, alpha=edge.color.alpha, curvature=edge.curve, arrow=arrow(length=unit(3,"pt"),type="closed")) 
	}else{
		gp <- gp + ggnetwork::geom_edges(color=edge.color, alpha=edge.color.alpha, curvature=edge.curve)
	}
	
	gp <- gp + ggnetwork::geom_nodes(aes(color=n.color,size=n.size), shape=node.shape)
	
	if(is.null(node.color.title)){
		node.color.title <- 'Node color'
	}
	if(is.null(zlim)){
		zlim <- range(df$n.color)
	}
	if(zlim[1] != zlim[2]){
		gp <- gp + scale_colour_gradientn(colors=xColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.5,nbin=5,draw.ulim=FALSE,draw.llim=FALSE))
	}else{
		gp <- gp + scale_colour_gradientn(colors=xColormap(colormap)(ncolors), guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.5))
	}
	
	if(is.null(node.size.title)){
		node.size.title <- 'Node size'
	}
	if(is.null(slim)){
		slim <- range(df$n.size)
	}
	if(slim[1] != slim[2]){
		gp <- gp + scale_size_continuous(limits=slim, range=node.size.range, guide=guide_legend(node.size.title,title.position="top",ncol=1))
	}else{
		gp <- gp + scale_size_continuous(limits=slim, range=node.size.range, guide=guide_legend(node.size.title,title.position="top",ncol=1))
	}
	
	gp <- gp + ggnetwork::theme_blank()
	gp <- gp + theme(text=element_text(family="sans")) + labs(title=title) + theme(plot.title=element_text(hjust=0.5), plot.margin=unit(rep(0,4),rep("lines",4)))
	
	if((zlim[1]!=zlim[2]) & (slim[1]!=slim[2])){
		gp <- gp + theme(legend.position="right")
    }else{
		if(slim[1]==slim[2]){
			gp <- gp + guides(size="none")
		}
		if(zlim[1]==zlim[2]){
			gp <- gp + guides(color="none")
		}
    }
    gp <- gp + theme(legend.title=element_text(size=8,face="bold"),legend.text=element_text(size=6))
    
    if(length(ls_ig)>1){
    	trait <- NULL
    	gp <- gp + facet_wrap(~trait)
		gp <- gp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=12,face="bold"), strip.placement="inside", panel.spacing=unit(0,"lines"))
    }
    
	if(sum(df$n.label=='')!=nrow(df)){
	
		###########
		StatNodes <- ggplot2::ggproto("StatNodes", ggplot2::Stat,
					compute_layer = function(data, scales, params) {
						if(all(c("xend", "yend") %in% names(data))) {
							unique(subset(data, select = c(-xend, -yend)))
						} else {
							unique(data)
						}
					}
				)
		## my_geom_nodetext_repel
		my_geom_nodetext_repel <- function (mapping = NULL, data = NULL, parse = FALSE, ..., box.padding = unit(0.25, "lines"), point.padding = unit(1e-06, "lines"), segment.size = 0.5, arrow = NULL, force = 1, max.iter = 2000, nudge_x = 0, nudge_y = 0, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE){
			ggplot2::layer(data = data, mapping = mapping, stat = StatNodes, 
				geom = ggrepel::GeomTextRepel, position = "identity", 
				show.legend = show.legend, inherit.aes = inherit.aes, 
				params = list(parse = parse, na.rm = na.rm, box.padding = box.padding, 
					point.padding = point.padding, 
					segment.size = segment.size, arrow = arrow, force = force, 
					max.iter = max.iter, nudge_x = nudge_x, nudge_y = nudge_y, 
					...)
			)
		}
		###########	
	
		if(length(unique(df$n.label.size))>1){
			#
			gp <- gp + my_geom_nodetext_repel(aes(label=n.label,size=n.label.size), lineheight=0.8, color=node.label.color, fontface="bold", alpha=node.label.alpha, box.padding=unit(0.5,"lines"), point.padding=unit(node.label.padding,"lines"), segment.alpha=0.5, segment.size=0.5, arrow=arrow(length=unit(node.label.arrow,'npc')), force=node.label.force)
		}else{
			if(unique(df$n.label.size)!=0){
				gp <- gp + my_geom_nodetext_repel(aes(label=n.label), lineheight=0.8, size=unique(df$n.label.size), color=node.label.color, fontface="bold", alpha=node.label.alpha, box.padding=unit(0.5,"lines"), point.padding=unit(node.label.padding,"lines"), segment.alpha=0.5, segment.size=0.5, arrow=arrow(length=unit(node.label.arrow,'npc')), force=node.label.force)
			}
		}
		
	}
    
    ####################
    ## append data_nodes
    if(1){
		df <- gp$data
		
		na.y <- NULL
		
		df_sub <- subset(df, is.na(na.y))
		ind <- match(colnames(df_sub), c('xend','yend','na.x','na.y'))
		df_sub <- df_sub[,is.na(ind)]
		df_sub <- df_sub[!duplicated(df_sub),]
		gp$data_nodes <- df_sub
    }
    
    invisible(gp)
}
