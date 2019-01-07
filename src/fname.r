library('tools')
fname = function(path) {
	#name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(path))
	name = file_path_sans_ext(path)
	ext  = file_ext(path)
	return(c(name,ext))
}