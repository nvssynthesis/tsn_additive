/*
  ==============================================================================

    logging.h
    Created: 14 Mar 2023 1:59:16pm
    Author:  Nicholas Solem

  ==============================================================================
*/

#ifndef NVS_LOGGING
#define NVS_LOGGING

#include "/Users/nicholassolem/development/thirdparty/Signalsmith/plot.h"
namespace nvs{
	class plotter{
	public:
		signalsmith::plot::PlotStyle customStyle() {
			signalsmith::plot::PlotStyle style;
			style.lineWidth = 1;
			style.valueSize = 9;
			style.fontAspectRatio = 1.1;
			style.fillOpacity = 0.6;
			style.tickH = style.tickV = 0;
			
			// Swap the first two colours, the second two dashes, and the 1st/3rd hatches
			std::swap(style.colours[0], style.colours[1]);
			std::swap(style.dashes[1], style.dashes[2]);
			std::swap(style.hatches[0], style.hatches[2]);
			
			style.cssSuffix = R"CSS(
				.svg-plot-value, .svg-plot-label {
					font-family: Verdana,sans-serif;
				}
				.svg-plot-axis {
					fill: #EEE;
				}
				.svg-plot-tick {
					stroke: #666;
					stroke-width: 0.75px;
				}
				.svg-plot-value {
					fill: #666;
					opacity: 0.8;
					font-weight: bold;
				}
				.svg-plot-major {
					stroke: #FFF;
					stroke-width: 1.5px;
				}
				.svg-plot-minor {
					stroke: #FFF;
					stroke-width: 0.75px;
					stroke-dasharray: none;
				}
			)CSS";
			// Minified version of `examples/wiggle.js`
			// This JS won't run inside an <img> tag - view the image itself, or embed it as <object>.
			style.scriptSrc = R"JS(var q=document,x=setTimeout,a=Math.random,r=q.querySelector("style");r.textContent='@import "/style/article/dist.css";'+r.textContent+'.svg-plot-value,.svg-plot-label{font-family:"Geraint Dense",Arial,sans-serif}';var t=[];q.querySelectorAll("path").forEach(function(e){var y=e.getAttribute("d");t.push(function(){var f=40*a(),u=!0,g=2*(a()-.5),h=2*(a()-.5),k=2*(a()-.5),l=2*(a()-.5),m,n,B=y.replace(/([0-9\.]+) ([0-9\.]+)/g,function(C,z,A){function p(b,c){b=parseFloat(b);c=parseFloat(c);if(!u){var d=b-m,v=c-n;d=Math.sqrt(d*
		d+v*v);if(20<d)return p(.5*(m+b),.5*(n+c)),p(b,c);f+=d;40<f&&(f=0,g=h,k=l,h=2*(a()-.5),l=2*(a()-.5))}u=!1;m=b;n=c;d=f/40;b+=g+(h-g)*d;c+=k+(l-k)*d;result+=" "+b+" "+c}result="";p(z,A);return result});e.setAttribute("d",B)})});var w=function(){t.forEach(function(e){e()});x(w,240*(.9+.2*a()))};w()
		)JS";
			return style;
		}
		template <typename T = float>
		void plotVector(std::vector<T> const &vec, std::string name = "default", const unsigned int maxTimesToPlot = 1) {
			if (timesPlotted < maxTimesToPlot)
			{
				name += std::to_string(timesPlotted);
				
				signalsmith::plot::Plot2D plot;

				// Customise the axes
				for (size_t i = 0; i < vec.size(); i += (vec.size() / 2)){
					plot.x.major(0).tick(i).label("time");
				}
				plot.y.major(0).minors(-1, 1).label("signal");

				// Add some data with `.add(x, y)`
				auto &figure = plot.line();
				
				for (size_t i = 0; i < vec.size(); ++i) {
					figure.add(i, vec[i]);
				}
				figure.label(name);

				plot.write(name + ".svg");
				timesPlotted++;
			}
		}
		template <typename T = float, size_t N = 8192>
		void plotArray(std::array<T, N>  &arr, std::string name = "default", const unsigned int maxTimesToPlot = 1, bool normalize = true) {
			if (timesPlotted < maxTimesToPlot)
			{
				name += std::to_string(timesPlotted);
				
				signalsmith::plot::Figure figure;
				auto &plot = figure.plot(1200, 400);

				// Customise the axes
				for (size_t i = 0; i < arr.size(); i += 32){
					plot.x.major(0).tick(i).label("time");
				}
				plot.y.major(0).minors(-1, 1).label("signal");

				// Add some data with `.add(x, y)`
				auto &thing = plot.line();
				
				if (normalize){
					typename std::array<T,N>::iterator result = std::max_element(arr.begin(), arr.end(), [](int a, int b) {
						return std::abs(a) < std::abs(b);
					});
					T norm = 1.0 / std::abs(*result);
					std::transform(arr.begin(), arr.end(), arr.begin(), [norm](auto i){
						return i * norm;
					});
				}
				
				for (size_t i = 0; i < arr.size(); ++i) {
					thing.add(i, arr[i]);
				}
				thing.label(name);

				figure.style = customStyle();
				figure.write(name + ".svg");
				timesPlotted++;
			}
		}
		void resetTimesPlotted(){
			timesPlotted = 0;
		}
	private:
		unsigned int timesPlotted {0};
	};
}

//void nvs::plotter::plotVector<float>(std::vector<float> vec, std::string name = "default", unsigned int maxTimesToPlot = 1);
#endif
