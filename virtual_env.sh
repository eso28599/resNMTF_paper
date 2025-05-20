if [ -f ".resnmtf_venv/bin/activate" ]; then
  source .resnmtf_venv/bin/activate
else
  # create virtual environment 
  if ! command -v pyenv &> /dev/null; then
    echo "pyenv is not installed. Please install pyenv first."
    exit 1
  fi

  if pyenv versions --bare | grep -qx 3.9.6; then
    echo "Python 3.9.6 is already installed via pyenv."
  else
    echo "Installing Python 3.9.6 via pyenv..."
    pyenv install 3.9.6 || {
      echo "Failed to install Python 3.9.6."
      exit 1
    }
    echo "Python 3.9.6 successfully installed."
  fi
  export PYENV_ROOT="$HOME/.pyenv"
  [[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
  eval "$(pyenv init - bash)"
  source ~/.bash_profile
  pyenv shell 3.9.6
  python -m venv .resnmtf_venv
  source .resnmtf_venv/bin/activate

  brew install openblas
  export PATH="/opt/homebrew/opt/openblas/bin:$PATH"
  export LDFLAGS="-L/opt/homebrew/opt/openblas/lib"
  export CPPFLAGS="-I/opt/homebrew/opt/openblas/include"

  python -m pip install -r requirements.txt
fi