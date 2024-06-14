import ternary as ternary

fig, tax = ternary.figure(scale=100)

fontsize = 15
offset = 0.05
tax.clear_matplotlib_ticks()
tax.gridlines(color="black", multiple=10)
tax.ticks(axis='lbr', linewidth=1, multiple=10)
tax.bottom_axis_label(r'a', fontsize=fontsize, color='blue')
tax.left_axis_label(r'b', fontsize=fontsize, color='red')
tax.right_axis_label(r'c', fontsize=fontsize, color='green')

tax.line((0,50,50), (50,0,50), color='red')
tax.line((30, 70, 0), (30, 0, 70), color='blue')
tax.line((30,20,50), (30,20,50), marker='s')
tax.line((0, 20, 80), (80, 20, 0), color='green')
# points are ordered in (a, c, b)


tax.get_axes().axis('off')
tax.show()

fig, tax = ternary.figure(scale=100)

fontsize = 15
offset = 0.05
tax.clear_matplotlib_ticks()
tax.gridlines(color="black", multiple=10)
tax.ticks(axis='lbr', linewidth=1, multiple=10)
tax.bottom_axis_label(r'a', fontsize=fontsize, color='blue')
tax.left_axis_label(r'b', fontsize=fontsize, color='red')
tax.right_axis_label(r'c', fontsize=fontsize, color='green')

# tax.line((0,50,50), (50,0,50), color='red')
# tax.line((30, 70, 0), (30, 0, 70), color='blue')
# tax.line((30,20,50), (30,20,50), marker='s')
# tax.line((0, 20, 80), (80, 20, 0), color='green')
# points are ordered in (a, c, b)


tax.get_axes().axis('off')
tax.show()