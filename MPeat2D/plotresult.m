%Plotting Peat Physical Properties
[xi, yi] = meshgrid( min(x_plot):1: max(x_plot),  min(y_plot):0.001: max(y_plot));
zi = griddata(x_plot,y_plot, (variable_plot), xi,yi);
contourf(xi,yi,(zi),100,'ShowText','off','edgecolor','none');
c = colorbar;
xlim([10 500]);
ylim([0 3]);
xlabel ('{\textbf{Horinzontal distance}} \boldmath$(m)$ ','Interpreter', 'Latex' );
ylabel ('{\textbf{Height}} \boldmath$(m)$ ','Interpreter', 'Latex' );
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontWeight', 'bold');