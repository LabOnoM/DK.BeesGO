import { QuartzComponent, QuartzComponentConstructor } from "./types"
// @ts-ignore
import script from "./scripts/followit.inline"

const FollowIt: QuartzComponent = () => {
  return (
    <>
      <div class="subscribe" style="text-align: left; margin: 0px 0;">
        <button id="open-subscribe-popup">
          <a class="button button--primary button--rounded button--xs" type="application/rss+xml">
            Subscribe by Email
          </a>
        </button>
      </div>
      <div
        id="subscribe-popup"
        style="display:none; position:fixed; z-index:9999; top:0; left:0; width:100vw; height:100vh; background:rgba(0,0,0,0.4);"
      >
        <div style="background:#fff; max-width:370px; width:90vw; margin:7vh auto 0 auto; border-radius:14px; box-shadow:0 10px 50px #3334; padding:32px 24px 16px 24px; position:relative;">
          <button
            id="close-subscribe-popup"
            style="position:absolute;top:10px;right:10px;font-size:22px;background:none;border:none;cursor:pointer;"
          >
            &times;
          </button>
          <form
            action="https://api.follow.it/subscription-form/bnRETXZPcTNXenFUU0lKWTcxTmNqeGhuZnY0cWQ0ZXI3REhFMkpObFZZVmJkdG01cmx4aGh1TlBtcDNNSkZsRFVTMWhMODFvRVcrRmdTbHJUVlh4S29SOERoT2FyTGdlNDJnMWE3R21nV3BYb290Rm1kdDhsZEFsbkJNRkFkczJ8MGI1NVJZdUhlSlpOazA3Y1pqRXlUQUYwSVhnQ0JXa21uQVEwRUdQQzByST0=/8"
            method="post"
          >
            <h5 style="color: #ad65c1; font-family: Arial; font-weight: bold; text-align: center;">
              Get new posts from DK.BeesGO by email:
            </h5>
            <input
              type="email"
              name="email"
              required
              placeholder="Enter your email"
              style="font-family: Arial; font-weight: normal;width:100%;height:40px; color: rgb(0, 0, 0); font-size: 14px; text-align: center; background-color: rgb(255, 255, 255); text-transform: none !important;"
            />
            <button
              type="submit"
              style="margin-top:10px;width:100%;height:40px;border-radius:6px;border:none;background:#6c78d5;color:#fff;font-weight:bold;font-size:16px;cursor:pointer;"
            >
              Subscribe
            </button>
          </form>
          <div style="text-align:center;font-size:12px;margin-top:8px;">
            Powered by{" "}
            <a href="https://follow.it" target="_blank" style="color:#4078c0;">
              follow.it
            </a>
          </div>
        </div>
      </div>
    </>
  )
}

FollowIt.afterDOMLoaded = script

export default (() => FollowIt) satisfies QuartzComponentConstructor
