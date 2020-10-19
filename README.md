# Description
Haystack was developed to assist gene hunting in our in-house variant database Sciobase&copy;. We presented Haystack as a poster during ISMB/ECCB '19. You can access our poster at [F1000](https://doi.org/10.7490/f1000research.1117489.1).

# Running Haystack
Haystack requires `R >= 3.6.x` to be run. It was implemented as a standalone Shiny Application and installs all dependencies on initial execution. To run Haystack, the Shiny application can be started with the Haystack source directory as working directory, for example using:

```
shiny:::runApp(port = 9999, launch.browser = FALSE)
```

You can use the file `example-variants.xlsx` from this repository to test Haystack's features.

# Contact

If you have any questions about Haystack, don't hesitate to reach out (mar.w@wwu.de)!
