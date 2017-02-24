import argparse
import ariba

def run(options):
    plotter = ariba.mic_plotter.MicPlotter(
      options.antibiotic,
      options.mic_file,
      options.summary_file,
      options.outprefix,
      use_hets=options.use_hets,
      main_title=options.main_title,
      plot_height=options.plot_height,
      plot_width=options.plot_width,
      log_y=not options.no_log_y,
      plot_types=options.plot_types,
      jitter_width=options.jitter_width,
      jitter_height=options.jitter_height,
      no_combinations=options.no_combinations,
      mic_values=options.mic_values,
      hlines=options.hlines,
      point_size=options.point_size,
      point_range=options.point_range,
      point_break=options.point_break,
      point_legend_x=options.point_legend_x,
      point_legend_y=options.point_legend_y,
      dot_size=options.dot_size,
      dot_outline=options.dot_outline,
      dot_y_text_size=options.dot_y_text_size,
      panel_heights=options.panel_heights,
      palette=options.palette,
      number_of_colours=options.number_of_colours
    )

    plotter.run()
