from meta_analysis import t_normalize, ternary_distance, ternary_subtract, ternary_dotproduct
import ternary
import math

if __name__ == "__main__":
    # ternary_points = [(30,40,30),(20,70,10),(2,90,3)] #,(15,80,5) #RH
    ternary_points = [(80, 10, 10), (10, 10, 80), (2, 90, 3)]  #LH

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
    main_line_slope = (ternary_points[0][1]-ternary_points[-1][1])/(ternary_points[0][0]-ternary_points[-1][0])

    main_line_eqn_output = main_line_slope*(max_point[0]-ternary_points[0][0]) - ternary_points[0][1] #point slope
    if main_line_slope < 0 and max_point[1] < main_line_eqn_output:
        max_distance = max_distance*-1
    elif main_line_slope > 0 and max_point[1] > main_line_eqn_output:
        max_distance = max_distance*-1
    # now max distance is the "curl" with sign correction
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