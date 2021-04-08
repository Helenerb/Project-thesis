library(wesanderson)
grand.budapest <- wes_palette("GrandBudapest1", n = 4)
grand.budapest.ligt <- c('#F8DFC0', '#FEB0B1', '#982B28', '#E39F76')

fantastic.fox <- wes_palette("GrandBudapest1", n=5)


data(iris)
head(iris, 6)

ggplot(data=iris) + geom_point(aes(x = Petal.Length, y = Sepal.Width, color = "width")) + 
  geom_point(aes(x = Petal.Length, y = Sepal.Length, color = "length")) + 
  scale_color_manual(name = "Type",
                     breaks = c("width", "length"),
                     values = c("width" = grand.budapest[1], "length" = grand.budapest[2]) )
