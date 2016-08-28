// Copyright ©2016 The gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package plotter

import (
	"fmt"
	"image/color"
	"math"
	"sort"

	"github.com/biogo/graphics/bezier"
	"github.com/gonum/plot"
	"github.com/gonum/plot/vg"
	"github.com/gonum/plot/vg/draw"
)

// A Sankey diagram presents stock and flow data as rectangles representing
// the amount of each stock and lines between the stocks representing the
// amount of each flow.
type Sankey struct {
	// Color specifies the default fill
	// colors for the stocks and flows. Colors can be
	// modified for individual stocks and flows.
	color.Color

	// StockBarWidth is the widths of the bars representing
	// the stocks. The default value is 15% larger than the
	// height of the stock label text.
	StockBarWidth vg.Length

	// LineStyle specifies the default border
	// line style for the stocks and flows. Styles can be
	// modified for individual stocks and flows.
	draw.LineStyle

	// TextStyle specifies the default stock label
	// text style. Styles can be modified for
	// individual stocks.
	draw.TextStyle

	flows []Flow

	// FlowStyle is a function that specifies the
	// background colar and border line style of the
	// flow based on its group name. The default
	// function uses the default Color and LineStyle
	// specified above for all groups.
	FlowStyle func(group string) (color.Color, draw.LineStyle)

	// stocks arranges the stocks by category.
	// The first key is the category and the seond
	// key is the label.
	stocks map[int]map[string]*stock
}

// stock represents the amount of a stock and its plotting order.
type stock struct {
	receptorValue, sourceValue float64
	label                      string
	category                   int
	order                      int

	// min represents the beginning of the plotting location
	// on the value axis.
	min float64

	// max is min plus the larger of receptorValue and sourceValue.
	max float64

	// sourceFlowPlaceholder and receptorFlowPlaceholder track
	// the current plotting location during
	// the plotting process.
	sourceFlowPlaceholder, receptorFlowPlaceholder float64
}

// A Flow represents the amount of an entity flowing between two stocks.
type Flow struct {
	// SourceLabel and ReceptorLabel are the labels
	// of the stocks that originate and receive the flow,
	// respectively.
	SourceLabel, ReceptorLabel string

	// SourceCategory and ReceptorCategory define
	// the locations on the category axis of the stocks that
	// originate and receive the flow, respectively. The
	// SourceCategory must be a lower number than
	// the ReceptorCategory.
	SourceCategory, ReceptorCategory int

	// Value represents the magnitute of the flow.
	// It must be greater than or equal to zero.
	Value float64

	// Group specifies the group that a flow belongs
	// to. It is used in assigning styles to groups
	// and creating legends.
	Group string
}

// NewSankey creates a new Sankey diagram with the specified
// flows and stocks.
func NewSankey(flows ...Flow) (*Sankey, error) {
	s := new(Sankey)

	s.stocks = make(map[int]map[string]*stock)

	s.flows = flows
	for i, f := range flows {
		// Here we make sure the stock categories are in the proper order.
		if f.SourceCategory >= f.ReceptorCategory {
			return nil, fmt.Errorf("plotter.NewSankey: Flow %d SourceCategory (%d) "+
				">= ReceptorCategory (%d)", i, f.SourceCategory, f.ReceptorCategory)
		}
		if f.Value < 0 {
			return nil, fmt.Errorf("plotter.NewSankey: Flow %d value (%g) < 0", i, f.Value)
		}

		// Here we initialize the stock holders.
		if _, ok := s.stocks[f.SourceCategory]; !ok {
			s.stocks[f.SourceCategory] = make(map[string]*stock)
		}
		if _, ok := s.stocks[f.ReceptorCategory]; !ok {
			s.stocks[f.ReceptorCategory] = make(map[string]*stock)
		}

		// Here we figure out the plotting order of the stocks.
		if _, ok := s.stocks[f.SourceCategory][f.SourceLabel]; !ok {
			s.stocks[f.SourceCategory][f.SourceLabel] = &stock{
				order:    len(s.stocks[f.SourceCategory]),
				label:    f.SourceLabel,
				category: f.SourceCategory,
			}
		}
		if _, ok := s.stocks[f.ReceptorCategory][f.ReceptorLabel]; !ok {
			s.stocks[f.ReceptorCategory][f.ReceptorLabel] = &stock{
				order:    len(s.stocks[f.ReceptorCategory]),
				label:    f.ReceptorLabel,
				category: f.ReceptorCategory,
			}
		}

		// Here we add the current value to the total value of the stocks
		s.stocks[f.SourceCategory][f.SourceLabel].sourceValue += f.Value
		s.stocks[f.ReceptorCategory][f.ReceptorLabel].receptorValue += f.Value
	}

	s.LineStyle = DefaultLineStyle
	s.LineStyle.Color = color.NRGBA{R: 0, G: 0, B: 0, A: 150}
	s.Color = color.NRGBA{R: 0, G: 0, B: 0, A: 100}

	fnt, err := vg.MakeFont(DefaultFont, DefaultFontSize)
	if err != nil {
		return nil, err
	}
	s.TextStyle = draw.TextStyle{
		Font:     fnt,
		Rotation: math.Pi / 2,
		XAlign:   draw.XCenter,
		YAlign:   draw.YCenter,
	}
	s.StockBarWidth = s.TextStyle.Font.Extents().Height * 1.15

	s.FlowStyle = func(_ string) (color.Color, draw.LineStyle) {
		return s.Color, s.LineStyle
	}

	return s, nil
}

// Plot implements the plot.Plotter interface.
func (s *Sankey) Plot(c draw.Canvas, plt *plot.Plot) {
	stocks := s.stockList()
	s.setStockMinMax(&stocks)

	trCat, trVal := plt.Transforms(&c)

	// Here we draw the flows.
	for _, f := range s.flows {
		startStock := s.stocks[f.SourceCategory][f.SourceLabel]
		endStock := s.stocks[f.ReceptorCategory][f.ReceptorLabel]
		catStart := trCat(float64(f.SourceCategory)) + s.StockBarWidth/2
		catEnd := trCat(float64(f.ReceptorCategory)) - s.StockBarWidth/2
		valStartLow := trVal(startStock.min + startStock.sourceFlowPlaceholder)
		valEndLow := trVal(endStock.min + endStock.receptorFlowPlaceholder)
		valStartHigh := trVal(startStock.min + startStock.sourceFlowPlaceholder + f.Value)
		valEndHigh := trVal(endStock.min + endStock.receptorFlowPlaceholder + f.Value)
		startStock.sourceFlowPlaceholder += f.Value
		endStock.receptorFlowPlaceholder += f.Value

		ptsLow := s.bezier(
			vg.Point{X: catStart, Y: valStartLow},
			vg.Point{X: catEnd, Y: valEndLow},
		)
		ptsHigh := s.bezier(
			vg.Point{X: catEnd, Y: valEndHigh},
			vg.Point{X: catStart, Y: valStartHigh},
		)

		color, lineStyle := s.FlowStyle(f.Group)

		// Here we fill the flow polygons.
		poly := c.ClipPolygonX(append(ptsLow, ptsHigh...))
		c.FillPolygon(color, poly)

		// Here we draw the flow edges.
		outline := c.ClipLinesX(ptsLow)
		c.StrokeLines(lineStyle, outline...)
		outline = c.ClipLinesX(ptsHigh)
		c.StrokeLines(lineStyle, outline...)
	}

	// Here we draw the stocks.
	for _, stk := range stocks {
		catLoc := trCat(float64(stk.category))
		if !c.ContainsX(catLoc) {
			continue
		}
		catMin, catMax := catLoc-s.StockBarWidth/2, catLoc+s.StockBarWidth/2
		valMin, valMax := trVal(stk.min), trVal(stk.max)

		// Here we fill the stock bars.
		pts := []vg.Point{
			{catMin, valMin},
			{catMin, valMax},
			{catMax, valMax},
			{catMax, valMin},
		}
		// poly := c.ClipPolygonX(pts) // This causes half of the bar to disappear. Is there a best practice here?
		c.FillPolygon(s.Color, pts) // poly)
		txtPt := vg.Point{X: (catMin + catMax) / 2, Y: (valMin + valMax) / 2}
		c.FillText(s.TextStyle, txtPt, stk.label)

		// Here we draw the bottom edge.
		pts = []vg.Point{
			{catMin, valMin},
			{catMax, valMin},
		}
		// outline := c.ClipLinesX(pts) // This causes half of the lines to disappear.
		c.StrokeLines(s.LineStyle, pts) //outline...)

		// Here we draw the top edge plus vertical edges where there are
		// no flows connected.
		pts = []vg.Point{
			{catMin, valMax},
			{catMax, valMax},
		}
		if stk.receptorValue < stk.sourceValue {
			y := trVal(stk.max - (stk.sourceValue - stk.receptorValue))
			pts = append([]vg.Point{{catMin, y}}, pts...)
		} else if stk.sourceValue < stk.receptorValue {
			y := trVal(stk.max - (stk.receptorValue - stk.sourceValue))
			pts = append(pts, vg.Point{X: catMax, Y: y})
		}
		//outline = c.ClipLinesX(pts)
		c.StrokeLines(s.LineStyle, pts) // outline...)
	}
}

// stockList returns a sorted list of the stocks in the diagram.
func (s *Sankey) stockList() []*stock {
	var stocks []*stock
	for _, ss := range s.stocks {
		for _, sss := range ss {
			stocks = append(stocks, sss)
		}
	}
	sort.Sort(stockSorter(stocks))
	return stocks
}

// stockSorter is a wrapper for a list of *stocks that implements
// sort.Interface.
type stockSorter []*stock

func (s stockSorter) Len() int      { return len(s) }
func (s stockSorter) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s stockSorter) Less(i, j int) bool {
	if s[i].category != s[j].category {
		return s[i].category < s[j].category
	}
	if s[i].order != s[j].order {
		return s[i].order < s[j].order
	}
	panic(fmt.Errorf("plotter: can't sort stocks:\n%+v\n%+v", s[i], s[j]))
}

// setStockMin sets the minimum and maximum values of the stock plotting locations.
func (s *Sankey) setStockMinMax(stocks *[]*stock) {
	var cat int
	var min float64
	for _, stk := range *stocks {
		stk.sourceFlowPlaceholder = 0
		stk.receptorFlowPlaceholder = 0
		if stk.category != cat {
			min = 0
		}
		cat = stk.category
		stk.min = min
		if stk.sourceValue > stk.receptorValue {
			stk.max = stk.min + stk.sourceValue
		} else {
			stk.max = stk.min + stk.receptorValue
		}
		min = stk.max
	}
}

// bezier creates a bezier curve between the begin and end points.
func (s *Sankey) bezier(begin, end vg.Point) []vg.Point {
	// directionOffsetFrac is the fraction of the distance between begin.X and
	// end.X for the bezier control points.
	const directionOffsetFrac = 0.3
	inPts := []vg.Point{
		begin,
		vg.Point{X: begin.X + (end.X-begin.X)*directionOffsetFrac, Y: begin.Y},
		vg.Point{X: begin.X + (end.X-begin.X)*(1-directionOffsetFrac), Y: end.Y},
		end,
	}
	curve := bezier.New(inPts...)

	// nPoints is the number of points for bezier interpolation.
	const nPoints = 20
	outPts := make([]vg.Point, nPoints)
	curve.Curve(outPts)
	return outPts
}

// DataRange implements the plot.DataRanger interface.
func (s *Sankey) DataRange() (xmin, xmax, ymin, ymax float64) {
	catMin := math.Inf(1)
	catMax := math.Inf(-1)
	for cat := range s.stocks {
		c := float64(cat)
		catMin = math.Min(catMin, c)
		catMax = math.Max(catMax, c)
	}

	valMin := math.Inf(1)
	valMax := math.Inf(-1)
	stocks := s.stockList()
	s.setStockMinMax(&stocks)
	for _, stk := range stocks {
		valMin = math.Min(valMin, stk.min)
		valMax = math.Max(valMax, stk.max)
	}
	return catMin, catMax, valMin, valMax
}

// GlyphBoxes implements the GlyphBoxer interface.
func (s *Sankey) GlyphBoxes(plt *plot.Plot) []plot.GlyphBox {
	stocks := s.stockList()
	s.setStockMinMax(&stocks)

	boxes := make([]plot.GlyphBox, 0, len(s.flows)+len(stocks))

	for _, stk := range stocks {
		b := plot.GlyphBox{
			X: plt.X.Norm(float64(stk.category)),
			Y: plt.Y.Norm((stk.min + stk.max) / 2),
			Rectangle: vg.Rectangle{
				Min: vg.Point{X: -s.StockBarWidth / 2},
				Max: vg.Point{X: s.StockBarWidth / 2},
			},
		}
		boxes = append(boxes, b)
	}
	return boxes
}

// Thumbnailers creates a group of objects that can be used to
// add legend entries for the different flow groups in this
// diagram, as well as the flow group labels that correspond to them.
func (s *Sankey) Thumbnailers() (legendLabels []string, thumbnailers []plot.Thumbnailer) {
	type empty struct{}
	flowGroups := make(map[string]empty)
	for _, f := range s.flows {
		flowGroups[f.Group] = empty{}
	}
	legendLabels = make([]string, len(flowGroups))
	thumbnailers = make([]plot.Thumbnailer, len(flowGroups))
	i := 0
	for g := range flowGroups {
		legendLabels[i] = g
		i++
	}
	sort.Strings(legendLabels)

	for i, g := range legendLabels {
		var thmb sankeyFlowThumbnailer
		thmb.Color, thmb.LineStyle = s.FlowStyle(g)
		thumbnailers[i] = plot.Thumbnailer(thmb)
	}
	return
}

// sankeyFlowThumbnailer implements the Thumbnailer interface
// for Sankey flow groups.
type sankeyFlowThumbnailer struct {
	draw.LineStyle
	color.Color
}

// Thumbnail fulfills the plot.Thumbnailer interface.
func (t sankeyFlowThumbnailer) Thumbnail(c *draw.Canvas) {
	// Here we draw the fill.
	pts := []vg.Point{
		{c.Min.X, c.Min.Y},
		{c.Min.X, c.Max.Y},
		{c.Max.X, c.Max.Y},
		{c.Max.X, c.Min.Y},
	}
	poly := c.ClipPolygonY(pts)
	c.FillPolygon(t.Color, poly)

	// Here we draw the upper border.
	pts = []vg.Point{
		{c.Min.X, c.Max.Y},
		{c.Max.X, c.Max.Y},
	}
	outline := c.ClipLinesY(pts)
	c.StrokeLines(t.LineStyle, outline...)

	// Here we draw the lower border.
	pts = []vg.Point{
		{c.Min.X, c.Min.Y},
		{c.Max.X, c.Min.Y},
	}
	outline = c.ClipLinesY(pts)
	c.StrokeLines(t.LineStyle, outline...)
}
