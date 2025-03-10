#!/bin/bash
HOME="/home/yuanj"
DOTFILES="/home/yuanj/Projects/dotfiles"

# ~/
cp $DOTFILES/condarc $HOME/.condarc
cp $DOTFILES/bashrc $HOME/.bashrc
cp $DOTFILES/gitconfig $HOME/.gitconfig
cp $DOTFILES/Rprofile $HOME/.Rprofile
cp $DOTFILES/SciTEUser.properties $HOME/.SciTEUser.properties

# ~/.local/
cp $DOTFILES/config/fcitx5/rime/default.custom.yaml $HOME/.local/share/fcitx5/rime/ -f
cp $DOTFILES/config/fcitx5/rime/user.yaml $HOME/.local/share/fcitx5/rime/ -f
cp $DOTFILES/config/fcitx5/rime/installation.yaml $HOME/.local/share/fcitx5/rime/ -f

# ~/.config/
cp $DOTFILES/config/ $HOME/.config/autostart -r
cp $DOTFILES/config/settings.json $HOME/.config/Code/User/settings.json -rf
cp $DOTFILES/config/ $HOME/.config/flameshot -r
cp $DOTFILES/config/ $HOME/.config/cinnamon -r

# /etc/
cp $DOTFILES/conf/pacman.conf /etc/
cp $DOTFILES/conf/lightdm.conf /etc/lightdm/
cp $DOTFILES/conf/64-language-selector-prefer.conf /etc/fonts/conf.d/

# /usr/lib/systemd/system/
cp $DOTFILES/systemd/rclone.service /usr/lib/systemd/system/
