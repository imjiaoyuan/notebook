#!/bin/bash
HOME="/home/yuanj"
DOTFILES="/home/yuanj/Projects/dotfiles"

# ~/
cp $HOME/.condarc $DOTFILES/condarc
cp $HOME/.bashrc $DOTFILES/bashrc
cp $HOME/.gitconfig $DOTFILES/gitconfig
cp $HOME/.Rprofile $DOTFILES/Rprofile
cp $HOME/.SciTEUser.properties $DOTFILES/SciTEUser.properties

# ~/.local/
cp $HOME/.local/share/fcitx5/rime/default.custom.yaml $DOTFILES/config/fcitx5/rime/ -f
cp $HOME/.local/share/fcitx5/rime/user.yaml $DOTFILES/config/fcitx5/rime/ -f
cp $HOME/.local/share/fcitx5/rime/installation.yaml $DOTFILES/config/fcitx5/rime/ -f

# ~/.config/
cp $HOME/.config/autostart $DOTFILES/config/ -r
cp $HOME/.config/Code/User/settings.json $DOTFILES/config/ -r
cp $HOME/.config/flameshot $DOTFILES/config/ -r
# cp $HOME/.config/cinnamon $DOTFILES/config/ -r

# /etc/
cp /etc/pacman.conf $DOTFILES/conf/
cp /etc/lightdm/lightdm.conf $DOTFILES/conf/
cp /etc/fonts/conf.d/64-language-selector-prefer.conf $DOTFILES/conf/

# /usr/lib/systemd/system/
cp /usr/lib/systemd/system/rclone.service $DOTFILES/systemd

# pacman -Qqen > packages-official.txt 
pacman -Qqen > $DOTFILES/misc/packages-official.txt

# pacman -Qqem > packages-archlinuxcn.txt
pacman -Qqem > $DOTFILES/misc/packages-archlinuxcn.txt
