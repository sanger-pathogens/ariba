import argparse
import ariba

def run(options):
    plotter = ariba.mic_plotter.MicPlotter(
      options.antibiotic,
      options.mic_file,
      options.summary_file,
      options.outprefix,
      main_title=options.main_title,
      plot_height=options.plot_height,
      plot_width=options.plot_width,
      log_y=not options.no_log_y,
      plot_types=options.plot_types,
      jitter_width=options.jitter_width,
      jitter_height=options.jitter_height,
      no_combinations=options.no_combinations
    )

    plotter.run()
