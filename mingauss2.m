	function MinimizingFunction = mingauss2(Point,data)

    MinimizingFunction=sum((data(:,2)-gaussian(data(:,1),Point(1),Point(2),Point(3))).^2);
