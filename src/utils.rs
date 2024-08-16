pub fn get_spinner() -> indicatif::ProgressBar{
    let bar = indicatif::ProgressBar::new_spinner();
    bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {pos} {per_sec}")
            .unwrap()
            .progress_chars("##-"),
    );
    bar
}