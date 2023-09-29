from meta_analysis import t_normalize, ternary_distance, ternary_subtract, ternary_dotproduct



if __name__ == "__main__":
    ternary_points = [(50,40,10),(2,95,3),(5,80,15)]


    def sig(x):
        return 1 if x > 0 else -1 if x < 0 else 0
    numerator_sum = 0
    if len(ternary_points) > 2:
        for i in range(1, len(ternary_points)):
            # get delta
            sign_factor = sig(ternary_dotproduct(
                ternary_subtract(ternary_points[-1],ternary_points[0]),
                ternary_subtract(ternary_points[i],ternary_points[i-1])
            )/(ternary_distance(ternary_points[-1],ternary_points[0])*ternary_distance(ternary_points[i],ternary_points[i-1])))
            # print(sig(sign_factor))

            numerator_sum = numerator_sum + ternary_distance(ternary_points[i],ternary_points[i-1])*sign_factor
    print(numerator_sum/ternary_distance(ternary_points[-1],ternary_points[0]))