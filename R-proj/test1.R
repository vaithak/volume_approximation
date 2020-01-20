ps = spectra_plot(2, 8, 5000, 800, 5)
 points1 = ps[1:2,1:10000]
points3 = ps[1:2,10801:11600]
points2 = ps[1:2,10001:10800]
 
 cbp1 <- c("grey3", "red4")
 ggplot(data.frame(x = c(points3[1,],points1[1,]), body = c(rep("P",800),
         rep("boundary",10000)), y = c(points3[2,],points1[2,])) , aes(x=x, y=y)) +
         scale_color_manual(values =cbp1) +
        geom_point(shape=20, aes(color=body)) +labs(x =" ", y = " ")