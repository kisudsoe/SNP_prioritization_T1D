fname = function(path) {
	name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(path))
	return(name)
}