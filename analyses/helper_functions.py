def exp_analyze(metric, formula, df, reg_output=True):

    import numpy as np
    import pandas as pd
    import datetime as dt
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    mod = smf.ols(formula=formula,
        data=df).fit(cov_type = 'HC1')

    #Extract parameters from regression
    coef = mod.params['variant']
    se = mod.bse['variant']
    pvalue = mod.pvalues['variant']
    means = df.groupby(['variant'])[metric].mean()
    control_mean = float(means[means.index == 0])
    rel_effect = 100*coef/control_mean

    lb = mod.conf_int(alpha=0.05, cols=None).loc['variant'][0]
    ub = mod.conf_int(alpha=0.05, cols=None).loc['variant'][1]
    ci_95 = abs(coef - mod.conf_int(alpha=0.05, cols=None).loc['variant'][0])
    ci_90 = abs(coef - mod.conf_int(alpha=0.1, cols=None).loc['variant'][0])
    ci_80 = abs(coef - mod.conf_int(alpha=0.2, cols=None).loc['variant'][0])

    #Print
    print("Analysis for {} metric:".format(metric))
    print("\nThe treatment coefficient was estimated to " + \
                          "have a value of {} ".format(round(coef, 4)) + \
                          "({}% relative to control) ".format(round(rel_effect, 2)) + \
                          "with a standard error of {}".format(round(se,4)) + \
                          " and a p-value of {}.".format(round(pvalue, 4)))
    print("\nThe 95% C.I. of the treatment variable coefficient " + \
                          "estimate ranges from {} to {}.\n".format(round(lb, 4), round(ub, 4)))
    print("\nThe control mean is {}".format(float(round(control_mean,2))))

    if reg_output == True:
        print('\n')
        print(mod.summary())

# Return dictionary with all needed variables for plotting

# Add in clustered and weighted analyses
    # if clustered=True:
    #
    #     #Clustered Regression
    #     mod = smf.ols(formula=formula, data=df).\
    #       fit(cov_type='cluster', cov_kwds={'groups': cluster})
    #
    # elif weighted=True:
    #
    #     #Weighted OLS
    #     mod = smf.wls(formula=formula, data=df, weights = weight).\
    #         fit(cov_type = 'HC1')


def exp_plot(metric, formula, df, suc_dir='up'):

    import math
    from scipy import stats
    import numpy as np
    import pandas as pd
    import datetime as dt
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    import matplotlib.pyplot as plt
    import seaborn as sns

    mod = smf.ols(formula=formula,
        data=df).fit(cov_type = 'HC1')

    #Extract parameters from regression
    coef = round(mod.params['variant'],4)
    pvalue = round(mod.pvalues['variant'],4)
    means = df.groupby(['variant'])[metric].mean()
    control_mean = round(float(means[means.index == 0]),4)
    rel_effect = round(100*coef/control_mean,4)

    lb = round(mod.conf_int(alpha=0.05, cols=None).loc['variant'][0],4)
    ub = round(mod.conf_int(alpha=0.05, cols=None).loc['variant'][1],4)
    ci_95 = round(abs(coef - mod.conf_int(alpha=0.05, cols=None).loc['variant'][0]),4)
    ci_90 = round(abs(coef - mod.conf_int(alpha=0.1, cols=None).loc['variant'][0]),4)
    ci_80 = round(abs(coef - mod.conf_int(alpha=0.2, cols=None).loc['variant'][0]),4)

    #Make figure
    fig, ax = plt.subplots(figsize=(5, 5))

    #Set Colors
    if pvalue > .05:
        s_color = 'gray'
    elif (coef < 0 and suc_dir == 'up') or (coef > 0 and suc_dir == 'down'):
        s_color = 'red'
    else:
        s_color = 'green'

    s_alpha = .5

    #Graphs
    ax.errorbar(y=0, x=coef, color=s_color, alpha = s_alpha, xerr=ci_95, lw=5, zorder = 1)
    ax.errorbar(y=0, x=coef, color=s_color, alpha = s_alpha, xerr=ci_90, lw=10, zorder = 1)
    ax.errorbar(y=0, x=coef, color=s_color, alpha = s_alpha, xerr=ci_80, lw=20, zorder = 1)
    ax.scatter(y=0, marker='o', s=120, x=coef, zorder = 10, color='gray')

    #Bounds of graph
    x_lower_bound = min(0,lb) - ci_95
    x_upper_bound = max(0,ub) + ci_95
    width = abs(x_upper_bound - x_lower_bound)

    y_lower_bound = -1
    y_upper_bound = 1
    height = abs(y_upper_bound - y_lower_bound)

    ax.set_xlim(x_lower_bound, x_upper_bound)
    ax.set_ylim(y_lower_bound, y_upper_bound)

    annotation_y = height/4 #y_lower_bound
    y_int = height*.03

    annotation_x = 1 #width*2.5 #x_lower_bound
    x_int = .42

    #Annotations
    plt.text(annotation_x, annotation_y + y_int, 'Abs Effect:', ha='left', fontsize = 12,transform=ax.transAxes)
    plt.text(annotation_x + x_int, annotation_y + y_int, coef, ha='right', fontsize = 12,transform=ax.transAxes)
    plt.text(annotation_x, annotation_y,'Rel Effect', ha='left', fontsize = 12,transform=ax.transAxes)
    plt.text(annotation_x + x_int, annotation_y,'{}%'.format(rel_effect),ha='right', fontsize = 12,transform=ax.transAxes)
    plt.text(annotation_x, annotation_y - y_int, 'p:', ha='left', fontsize = 12,transform=ax.transAxes)
    plt.text(annotation_x + x_int, annotation_y - y_int, pvalue, ha='right', fontsize = 12,transform=ax.transAxes)

    #title
    ax.set_title('Effect of treatment on {}\n'.format(metric), fontsize = 16)

    #baseline
    ax.axvline(x=0, linestyle='--', color='gray', linewidth=1, zorder = 0)

    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.set_yticklabels('')
    ax.set_frame_on(False)

    plt.show()
