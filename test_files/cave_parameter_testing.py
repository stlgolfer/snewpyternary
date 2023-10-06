from meta_analysis import t_normalize, ternary_distance, ternary_subtract, ternary_dotproduct
import ternary
import math

if __name__ == "__main__":
    ternary_points = [(30,40,30),(20,70,10),(15,80,5),(2,90,3)]

    figure, tax = ternary.figure(scale=100)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="black", multiple=10)
    # data is organized in top, right, left

    # heatmap stuff
    # if not heatmap == None:
    #     tax.heatmap(heatmap)

    tax.scatter(points=ternary_points)
    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')  # disables regular matlab plot axe

    # now choose the point that is farthest away. should be O(n) so should be ok
    max_point = None
    max_distance = 0
    main_line = ternary_distance(ternary_points[0], ternary_points[-1])
    for point in ternary_points[1:-1]:
        a = ternary_distance(point, ternary_points[-1])
        b = ternary_distance(point, ternary_points[0])
        s = (a + main_line + b)/2 # Heron's semi-perimeter
        h = math.sqrt(4*s*(s-a)*(s-b)*(s-main_line)/main_line**2)
        if h > max_distance:
            max_distance = h
            max_point = point

    # this is probably a bad way of doing this, but let's measure the sign by first
    # find the axis that the main line is "most" on. For example, if the end points are
    # basically stacked on top of each other, see if the max point is to the right or left
    # since most of the models work this way, it should be ok. actually just use whether
    # it's right or left of the end point (so somewhere in the middle)
    if max_point[0] < ternary_points[-1][0]:
        max_distance = -1*max_distance
    # now max distance is the "curl"
    print(max_distance)
    tax.get_axes().text(0,-10,rf'$\iota={round(max_distance,3)}$')
    tax.show()

    #region old way of doing this
    # def sig(x):
    #     return 1 if x > 0 else -1 if x < 0 else 0
    # numerator_sum = 0
    # if len(ternary_points) > 2:
    #     for i in range(1, len(ternary_points)):
    #         # get delta
    #         sign_factor = sig(ternary_dotproduct(
    #             ternary_subtract(ternary_points[-1],ternary_points[0]),
    #             ternary_subtract(ternary_points[i],ternary_points[i-1])
    #         )/(ternary_distance(ternary_points[-1],ternary_points[0])*ternary_distance(ternary_points[i],ternary_points[i-1])))
    #         # print(sig(sign_factor))
    #
    #         numerator_sum = numerator_sum + ternary_distance(ternary_points[i],ternary_points[i-1])*sign_factor
    # print(numerator_sum/ternary_distance(ternary_points[-1],ternary_points[0]))
    #endregion