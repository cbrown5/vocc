# To do for vocc R package
2018-01-12

This list is for developers (Dave Schoeman and Chris Brown currently)

## 1

Fix zero gradient problem. Some areas, e.g. Strait of Gibraltar have zero spatial gradients resulting in NA velocity values. This is an issue for species trajectory models, because it creates artificial barriers.

## 2

There are some differences between new and old code in cul-de-sacs, due to gradients being slightly different. My new code tends to calculate much greater values of vocc in these places.
Overall this won't affect maps much, but would be good to know why difference occurs.
See "data-raw/compare_new_old_functions.R" for a comparison.

## 3
Create multivariate version of the Hamann function
