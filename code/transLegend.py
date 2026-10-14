from matplotlib.transforms import BboxTransformTo, Bbox

def transLegend(legend):
	fig = legend.get_figure()
	ax = legend.axes
	fig.canvas.draw()

	bbox_display = legend.get_window_extent()
	bbox_axesfrac = Bbox(ax.transAxes.inverted().transform(bbox_display))

	# unit square -> axes-fraction bbox -> display (dynamically, via ax.transAxes)
	return BboxTransformTo(bbox_axesfrac) + ax.transAxes

def get_transform_aspect(transform):
	p0 = transform.transform((0, 0))
	px = transform.transform((1, 0))
	py = transform.transform((0, 1))
	xscale = px[0] - p0[0]
	yscale = py[1] - p0[1]
	return xscale / yscale
