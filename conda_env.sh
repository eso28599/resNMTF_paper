# Check if 'resnmtf_env' is already activated
# export PATH=/opt/anaconda3/bin:$PATH
if [[ "$CONDA_DEFAULT_ENV" == "resnmtf_env" ]]; then
    echo "Environment 'resnmtf_env' is already activated."
elif conda env list | grep -q "resnmtf_env"; then
    # Environment exists but is not activated
    echo "Environment 'resnmtf_env' exists. Activating it..."
    conda init zsh
    source ~/.zshrc
    conda activate resnmtf_env
else
    # Environment does not exist, so create and activate it
    echo "Environment 'resnmtf_env' does not exist. Creating and activating..."
    conda env create -f environment.yml
    conda activate resnmtf_env
fi